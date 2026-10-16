// Geant4 truth energy spectrum for the LaBr3:Ce crystal, split by primary source.
// The PhotonGunPlusLa138 simulation feeds two independent generators into one g4run:
//   1) PhotonGun  - the 1.809 MeV signal gamma, fired from the STM aperture (z=40404),
//                   well upstream of the crystal.
//   2) La138Source - intrinsic La-138 decay products emitted from uniform vertices
//                    INSIDE the crystal volume (z ~ 40891 +/- 38 mm).
// Both sets of primaries are transported by Geant4 and deposit energy in the same
// StepPointMC collection (g4run:STMDet), so the digitized waveform mixes them.
// To recover the truth spectrum of each contribution separately, this analyzer traces
// every energy-depositing step back to its PRIMARY and classifies it with
// LaBrTruth::classifyPrimary (LaBrTruthHelpers.hh): by generator tag (La-138, Ac-227), else
// inside/outside the crystal. "Gun" means any EXTERNAL source (photon gun or calibration
// source) wherever it is placed. Per event it sums the light-equivalent energy (alpha steps
// scaled by alphaQuenchingFactor, as in the digitizer) deposited inside the LaBr crystal into
// four branches (total, external, La-138, Ac-227) in TruthLaBr/ttree.
//
// It also writes TruthLaBr/decays: ONE ENTRY PER DECAY (or photon-gun photon) that deposited
// energy, with the decay's total light-equivalent deposit, its source and the time of its
// first deposit. Products of one decay are grouped by their common vertex and start time
// (LaBrTruth::decayKey), so a La-138 EC decay's gamma, X-ray and Auger electron are one entry.
// This is the pure Geant4 truth spectrum per decay; it does not depend on the pileup chain. The selection matches the
// digitizer's: the LaBr-side x cut (to reject the HPGe crystal sharing the collection) plus
// the LaBr crystal-cylinder bounds. No time cut is applied. This keeps the truth spectrum
// over the same crystal volume the digitizer sees, so the two are directly comparable.
// Original author: Sena Clarke

// stdlib includes
#include <cmath>
#include <map>

// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// fhicl includes
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline includes
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/STMMC/inc/LaBrTruthHelpers.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

// ROOT includes
#include "art_root_io/TFileService.h"
#include "TTree.h"

namespace mu2e {
  class LaBrTruthSpectrum : public art::EDAnalyzer {
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config {
      fhicl::Atom<art::InputTag> StepPointMCsTag{Name("StepPointMCsTag"), Comment("StepPointMCs from g4run, e.g. g4run:STMDet")};
      fhicl::Atom<double> alphaQuenchingFactor{Name("alphaQuenchingFactor"), Comment("Alpha/electron light ratio at 6 MeV; must match the digitizer. Default 1."), 1.0};
      fhicl::Atom<double> alphaQuenchingSlope{Name("alphaQuenchingSlope"), Comment("Change of the alpha quenching factor per MeV [1/MeV]; must match the digitizer. Default 0."), 0.0};
      fhicl::Atom<bool> makeStepTree{Name("makeStepTree"), Comment("Write the per-step deposit TTree (x,y,z,edep,source) for the deposit-weighted spatial map. One entry per kept step, so it grows with statistics; off by default."), false};
    };
    using Parameters = art::EDAnalyzer::Table<Config>;
    explicit LaBrTruthSpectrum(const Parameters& conf);
    void analyze(const art::Event& event) override;
    void endJob() override;
  private:
    art::ProductToken<StepPointMCCollection> StepPointMCsToken;
    double alphaQF = 1.0, alphaQFSlope = 0.0;
    bool makeStepTree = false;

    // Accept LaBr-side deposits only (x > this), rejecting the HPGe crystal which shares
    // the g4run:STMDet collection. HPGe centre x ~ -3986, LaBr ~ -3863, midpoint ~ -3925.
    const double labrSideXcut = -3925.0;

    // LaBr crystal geometry [mm], identical to the LaBr digitizer's selection, so the truth
    // spectrum is restricted to the same crystal volume the digitizer sees.
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

    // Per-event TTree variables (energies in keV)
    TTree* ttree = nullptr;
    uint eventId = 0;
    double edepTotal = 0, edepGun = 0, edepLa138 = 0, edepAc227 = 0;

    // Per-step TTree: one entry per energy-depositing step kept by the selection above,
    // carrying its world position and ionizing deposit so a deposit-weighted x-y map can be
    // built (TH2D filled with stepEdep as the weight). stepSource classifies the primary:
    // 0 = external (gun/source), 1 = La-138, 2 = Ac-227. stepPrimaryId is the originating primary's SimParticle id,
    // so (eventId, stepPrimaryId) groups together every deposit descended from one primary
    // photon - allowing a single photon's interaction history to be drawn on its own.
    TTree* stepTree = nullptr;

    // Per-decay TTree
    TTree* decayTree = nullptr;
    double decayEnergy = 0, decayTime = 0;
    int decaySource = 0;
    struct DecayAccum { double energy = 0; double tMin = 0; bool seen = false; };
    double stepX = 0, stepY = 0, stepZ = 0, stepEdep = 0;
    int stepSource = 0;
    uint stepPrimaryId = 0;

    // Diagnostic counters
    uint64_t n_events = 0, n_steps_gun = 0, n_steps_la138 = 0, n_steps_ac227 = 0;
  };

  LaBrTruthSpectrum::LaBrTruthSpectrum(const Parameters& conf) :
    art::EDAnalyzer(conf),
    StepPointMCsToken(consumes<StepPointMCCollection>(conf().StepPointMCsTag())),
    alphaQF(conf().alphaQuenchingFactor()),
    alphaQFSlope(conf().alphaQuenchingSlope()),
    makeStepTree(conf().makeStepTree()) {
      art::ServiceHandle<art::TFileService> tfs;
      ttree = tfs->make<TTree>("ttree", "LaBr Geant4 truth energy spectrum by source");
      ttree->Branch("eventId", &eventId, "eventId/i");
      ttree->Branch("edepTotal", &edepTotal, "edepTotal/D"); // keV, all sources
      ttree->Branch("edepGun", &edepGun, "edepGun/D");       // keV, external source only
      ttree->Branch("edepLa138", &edepLa138, "edepLa138/D"); // keV, La-138 only
      ttree->Branch("edepAc227", &edepAc227, "edepAc227/D"); // keV, Ac-227 chain only (light-equivalent)

      decayTree = tfs->make<TTree>("decays", "LaBr Geant4 truth deposit per decay");
      decayTree->Branch("eventId", &eventId, "eventId/i");
      decayTree->Branch("energy", &decayEnergy, "energy/D"); // keV, light-equivalent, whole decay
      decayTree->Branch("source", &decaySource, "source/I"); // 0 = external, 1 = La-138, 2 = Ac-227
      decayTree->Branch("time", &decayTime, "time/D");       // ns within the micropulse, first deposit

      if (makeStepTree) {
        stepTree = tfs->make<TTree>("steps", "LaBr Geant4 truth deposits, one entry per step");
        stepTree->Branch("eventId", &eventId, "eventId/i");
        stepTree->Branch("x", &stepX, "x/D");           // world x [mm]
        stepTree->Branch("y", &stepY, "y/D");           // world y [mm]
        stepTree->Branch("z", &stepZ, "z/D");           // world z [mm]
        stepTree->Branch("edep", &stepEdep, "edep/D");  // ionizing deposit [keV]
        stepTree->Branch("source", &stepSource, "source/I"); // 0 = external, 1 = La-138, 2 = Ac-227
        stepTree->Branch("primaryId", &stepPrimaryId, "primaryId/i"); // originating primary's SimParticle id
      }
    };

  void LaBrTruthSpectrum::analyze(const art::Event& event) {
    eventId = event.id().event();
    ++n_events;
    edepTotal = edepGun = edepLa138 = edepAc227 = 0;
    std::map<LaBrTruth::DecayKey, DecayAccum> byDecay;

    auto const& steps = event.getProduct(StepPointMCsToken);
    for (const StepPointMC& step : steps) {
      if (step.ionizingEdep() <= 0)
        continue;

      // Keep only LaBr-side deposits (reject the HPGe crystal), then require the step to
      // lie within the LaBr crystal cylinder, matching the digitizer's selection so the
      // truth spectrum covers the same volume.
      const CLHEP::Hep3Vector& world = step.position();
      if (world.x() < labrSideXcut)
        continue;
      double localX = world.x() - crystalCentreX;
      double localY = world.y() - crystalCentreY;
      double localZ = (world.z() - crystalCentreZ) + (crystalL / 2); // front face at z=0
      double localR = std::sqrt(localX * localX + localY * localY);
      if (localZ < -stepPositionTolerance || localZ > maxZ || localR > maxR)
        continue;

      // Trace this step to its primary (walking the parent chain explicitly:
      // SimParticle::originParticle() returns ordinary secondaries themselves) and classify it.
      const SimParticle* primary = LaBrTruth::findPrimary(step);
      if (!primary)
        continue;
      double edep_keV = LaBrTruth::lightEquivalentEdep(step, alphaQF, alphaQFSlope) * 1e3;
      edepTotal += edep_keV;
      const int source = LaBrTruth::classifyPrimary(*primary, crystal);
      DecayAccum& d = byDecay[LaBrTruth::decayKey(*primary, source)];
      d.energy += edep_keV;
      if (!d.seen || step.time() < d.tMin) d.tMin = step.time();
      d.seen = true;
      if (source == STMMCPulse::Ac227) {
        edepAc227 += edep_keV;
        ++n_steps_ac227;
      } else if (source == STMMCPulse::La138) {
        edepLa138 += edep_keV;
        ++n_steps_la138;
      } else {
        edepGun += edep_keV;
        ++n_steps_gun;
      }

      // Record this deposit for the deposit-weighted spatial map (opt-in).
      if (makeStepTree) {
        stepX = world.x();
        stepY = world.y();
        stepZ = world.z();
        stepEdep = edep_keV;
        stepSource = source;
        stepPrimaryId = primary->id().asUint();
        stepTree->Fill();
      }
    }

    ttree->Fill();

    for (const auto& kv : byDecay) {
      decaySource = std::get<0>(kv.first);
      decayEnergy = kv.second.energy;
      decayTime = kv.second.tMin;
      decayTree->Fill();
    }
  }

  void LaBrTruthSpectrum::endJob() {
    mf::LogInfo("LaBrTruthSpectrum")
      << "events: " << n_events
      << ", external-source steps: " << n_steps_gun
      << ", La-138 steps: " << n_steps_la138
      << ", Ac-227 steps: " << n_steps_ac227;
  }
}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::LaBrTruthSpectrum)
