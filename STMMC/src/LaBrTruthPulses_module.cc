// Stage 1 of the STM pileup-truth chain.
//
// Turns the Geant4 energy deposits in the LaBr crystal into one STMMCPulse per
// DECAY (or per photon-gun photon): the pulse is all the ionizing energy deposited by
// every primary of that decay, placed at the earliest of their deposit times (the
// leading edge the pulse finder triggers on). A decay emits several primaries (La-138
// EC: 1436 keV gamma + X-ray + Auger electron; Co-60: two gammas) that arrive within
// nanoseconds and form ONE pulse in the detector, so primaries are grouped by their
// common generator vertex and start time. Grouping per primary instead would split the
// 1473 keV sum peak into 1436 + 37 keV pulses, count particles rather than decays, and
// make TagPileupPulses tag a decay's own products as pileup of each other. Times are converted to an ABSOLUTE axis shared across micropulses
// (art events) by adding a per-event origin that advances by one micropulse each
// event, so a pulse at the end of event N and one at the start of event N+1 land
// on a single continuous timeline. That is what lets the downstream tagger
// (TagPileupPulses) detect pileup across micropulse boundaries.
//
// The crystal selection (LaBr-side x cut + crystal-cylinder bounds) and the
// primary tracing (walking the SimParticle parent chain) match
// LaBrTruthSpectrum_module.cc exactly, so the pulses cover the same volume the
// LaBr digitizer sees. No pileup decision is made here; isPileup is left false.
//
// Source classification and alpha light quenching come from LaBrTruthHelpers.hh:
// source 0 = external (photon gun / calibration source), 1 = La-138, 2 = Ac-227 chain.
// Pulse energies are light-equivalent (alpha steps scaled by alphaQuenchingFactor), so
// they match what the digitizer produces.
//
// Original author: Sena Clarke

// stdlib includes
#include <algorithm>
#include <cmath>
#include <map>

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// fhicl includes
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline includes
#include "Offline/MCDataProducts/inc/STMMCPulse.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/STMMC/inc/LaBrTruthHelpers.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {
  class LaBrTruthPulses : public art::EDProducer {
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config {
      fhicl::Atom<art::InputTag> StepPointMCsTag{Name("StepPointMCsTag"), Comment("StepPointMCs from g4run, e.g. g4run:STMDet")};
      fhicl::Atom<double> micropulseTime{Name("micropulseTime"), Comment("Micropulse (event) length [ns]; must match the digitizer's micropulseTime. Default 1695."), 1695.0};
      fhicl::Atom<double> alphaQuenchingFactor{Name("alphaQuenchingFactor"), Comment("Alpha/electron light ratio at 6 MeV; must match the digitizer. Default 1."), 1.0};
      fhicl::Atom<double> alphaQuenchingSlope{Name("alphaQuenchingSlope"), Comment("Change of the alpha quenching factor per MeV [1/MeV]; must match the digitizer. Default 0."), 0.0};
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit LaBrTruthPulses(const Parameters& conf);
    void produce(art::Event& event) override;
    void endJob() override;
  private:
    art::ProductToken<StepPointMCCollection> StepPointMCsToken;
    double micropulseTime = 1695.0;
    double alphaQF = 1.0, alphaQFSlope = 0.0;

    // Running origin of the current micropulse on the absolute time axis [ns].
    // Advanced by micropulseTime every event, whether or not it had a deposit,
    // so the timeline stays aligned with the sequential EmptyEvent micropulses.
    double _timeOrigin = 0.0;

    // LaBr crystal geometry [mm] and selection, identical to LaBrTruthSpectrum.
    const double labrSideXcut = -3925.0;
    const double crystalCentreX = -3863.3600;
    const double crystalCentreY = 0.0;
    const double crystalCentreZ = 40891.5000;
    const double crystalL = 76.2;
    const double crystalR = 38.1;
    const double stepPositionTolerance = 0.1;
    const double maxR = crystalR + stepPositionTolerance;
    const double maxZ = crystalL + stepPositionTolerance;
    const LaBrTruth::CrystalCylinder crystal{crystalCentreX, crystalCentreY, crystalCentreZ,
                                             crystalR, crystalL / 2};

    // Per-decay accumulator within one event.
    struct Accum {
      double energy = 0.0;   // summed ionizing deposit [keV]
      double tMin = 0.0;     // earliest deposit time [ns], local to the micropulse
      int source = STMMCPulse::unknown;
      unsigned primaryId = 0; // lowest SimParticle id among the decay's primaries
      bool seen = false;
    };

    // Diagnostics
    uint64_t n_events = 0, n_pulses = 0;
  };

  LaBrTruthPulses::LaBrTruthPulses(const Parameters& conf) :
    art::EDProducer(conf),
    StepPointMCsToken(consumes<StepPointMCCollection>(conf().StepPointMCsTag())),
    micropulseTime(conf().micropulseTime()),
    alphaQF(conf().alphaQuenchingFactor()),
    alphaQFSlope(conf().alphaQuenchingSlope()) {
      produces<STMMCPulseCollection>();
    }

  void LaBrTruthPulses::produce(art::Event& event) {
    ++n_events;
    const unsigned eventId = event.id().event();
    std::unique_ptr<STMMCPulseCollection> out(new STMMCPulseCollection);

    std::map<LaBrTruth::DecayKey, Accum> byDecay;

    auto const& steps = event.getProduct(StepPointMCsToken);
    for (const StepPointMC& step : steps) {
      if (step.ionizingEdep() <= 0)
        continue;

      // LaBr crystal selection (reject HPGe side, require in-cylinder), matching
      // the digitizer so pulses cover the same volume it digitizes.
      const CLHEP::Hep3Vector& world = step.position();
      if (world.x() < labrSideXcut)
        continue;
      double localX = world.x() - crystalCentreX;
      double localY = world.y() - crystalCentreY;
      double localZ = (world.z() - crystalCentreZ) + (crystalL / 2); // front face at z=0
      double localR = std::sqrt(localX * localX + localY * localY);
      if (localZ < -stepPositionTolerance || localZ > maxZ || localR > maxR)
        continue;

      const SimParticle* primary = LaBrTruth::findPrimary(step);
      if (!primary)
        continue;

      const unsigned pid = primary->id().asUint();
      const int source = LaBrTruth::classifyPrimary(*primary, crystal);
      Accum& a = byDecay[LaBrTruth::decayKey(*primary, source)];
      double edep_keV = LaBrTruth::lightEquivalentEdep(step, alphaQF, alphaQFSlope) * 1e3;
      if (!a.seen || step.time() < a.tMin)
        a.tMin = step.time();
      a.energy += edep_keV;
      a.source = source;
      a.primaryId = a.seen ? std::min(a.primaryId, pid) : pid;
      a.seen = true;
    }

    for (const auto& kv : byDecay) {
      const Accum& a = kv.second;
      double tAbs = _timeOrigin + a.tMin;
      out->emplace_back(tAbs, a.energy, a.primaryId, a.source, eventId);
      ++n_pulses;
    }

    _timeOrigin += micropulseTime;
    event.put(std::move(out));
  }

  void LaBrTruthPulses::endJob() {
    mf::LogInfo("LaBrTruthPulses")
      << "events: " << n_events << ", truth pulses emitted: " << n_pulses;
  }
}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::LaBrTruthPulses)
