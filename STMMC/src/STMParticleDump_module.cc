// Writes the particles crossing one or more virtual detectors to a flat TTree that
// STMParticleSource (mode "dataset") resamples as the input of the LaBr pipeline.
// Use it on the output of an upstream stage (e.g. BeamToVD / NeutralsToVD, tag
// "compressDetStepMCsSTM:virtualdetector"), so the photons, neutrons, electrons... that reach
// the STM carry their real spectra, directions, times and per-event correlations.
//
// Tree "particles", one entry per VD crossing:
//   eventId (L)             sequential input-event index, groups particles of one upstream event
//   run, subRun, event (I)  art ID of that event
//   virtualdetectorId (I), pdgId (I), trackId (L)
//   x, y, z [mm], px, py, pz [MeV/c], time [ns], E [MeV, kinetic]
// Tree "info", one entry: nEvents (L) = input events read, for normalisation:
//   particles per upstream event = particles->GetEntries() / nEvents.
// Original author: Sena Clarke

// stdlib includes
#include <cmath>
#include <map>
#include <set>
#include <utility>
#include <vector>

// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"

// exception handling
#include "cetlib_except/exception.h"

// fhicl includes
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline includes
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"

// ROOT includes
#include "art_root_io/TFileService.h"
#include "TTree.h"

namespace mu2e {
  class STMParticleDump : public art::EDAnalyzer {
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config {
      fhicl::Atom<art::InputTag> StepPointMCsTag{Name("StepPointMCsTag"), Comment("Virtual-detector StepPointMCs, e.g. compressDetStepMCsSTM:virtualdetector")};
      fhicl::Sequence<int> virtualDetectorIds{Name("virtualDetectorIds"), Comment("VD ids to keep, e.g. [101]. Empty = all."), std::vector<int>{}};
      fhicl::Sequence<int> pdgIds{Name("pdgIds"), Comment("PDG codes to keep. Empty = all."), std::vector<int>{}};
      fhicl::Atom<bool> firstCrossingOnly{Name("firstCrossingOnly"), Comment("Keep only the earliest crossing of each track through each VD (drops backscatter re-crossings)"), true};
      fhicl::Atom<bool> forwardOnly{Name("forwardOnly"), Comment("Keep only crossings with pz > 0 (moving downstream)"), false};
    };
    using Parameters = art::EDAnalyzer::Table<Config>;
    explicit STMParticleDump(const Parameters& conf);
    void analyze(const art::Event& event) override;
    void endJob() override;

  private:
    art::ProductToken<StepPointMCCollection> _stepsToken;
    std::set<int> _vds, _pdgs;
    bool _firstOnly, _forwardOnly;
    GlobalConstantsHandle<ParticleDataList> _pdt;

    TTree* _particles = nullptr;
    TTree* _info = nullptr;
    Long64_t _eventId = 0, _trackId = 0, _nEvents = 0;
    int _run = 0, _subRun = 0, _event = 0, _vdId = 0, _pdgId = 0;
    double _x = 0, _y = 0, _z = 0, _px = 0, _py = 0, _pz = 0, _time = 0, _E = 0;
    std::map<int, unsigned long> _pdgCounts;
  };

  STMParticleDump::STMParticleDump(const Parameters& conf) :
    art::EDAnalyzer(conf),
    _stepsToken(consumes<StepPointMCCollection>(conf().StepPointMCsTag())),
    _firstOnly(conf().firstCrossingOnly()),
    _forwardOnly(conf().forwardOnly()) {
      for (int id : conf().virtualDetectorIds()) _vds.insert(id);
      for (int id : conf().pdgIds()) _pdgs.insert(id);
      art::ServiceHandle<art::TFileService> tfs;
      _particles = tfs->make<TTree>("particles", "Particles crossing the virtual detectors");
      _particles->Branch("eventId", &_eventId, "eventId/L");
      _particles->Branch("run", &_run, "run/I");
      _particles->Branch("subRun", &_subRun, "subRun/I");
      _particles->Branch("event", &_event, "event/I");
      _particles->Branch("virtualdetectorId", &_vdId, "virtualdetectorId/I");
      _particles->Branch("pdgId", &_pdgId, "pdgId/I");
      _particles->Branch("trackId", &_trackId, "trackId/L");
      _particles->Branch("x", &_x, "x/D");
      _particles->Branch("y", &_y, "y/D");
      _particles->Branch("z", &_z, "z/D");
      _particles->Branch("px", &_px, "px/D");
      _particles->Branch("py", &_py, "py/D");
      _particles->Branch("pz", &_pz, "pz/D");
      _particles->Branch("time", &_time, "time/D");
      _particles->Branch("E", &_E, "E/D");
      _info = tfs->make<TTree>("info", "Normalisation");
      _info->Branch("nEvents", &_nEvents, "nEvents/L");
    };

  void STMParticleDump::analyze(const art::Event& event) {
    _run = event.run();
    _subRun = event.subRun();
    _event = event.event();
    auto const& steps = event.getProduct(_stepsToken);

    // Select, and with firstCrossingOnly keep the earliest step per (track, VD)
    std::map<std::pair<unsigned long, int>, const StepPointMC*> first;
    std::vector<const StepPointMC*> kept;
    for (const StepPointMC& step : steps) {
      const int vd = static_cast<int>(step.virtualDetectorId());
      if (!_vds.empty() && !_vds.count(vd)) continue;
      const int pdg = step.simParticle()->pdgId();
      if (!_pdgs.empty() && !_pdgs.count(pdg)) continue;
      if (_forwardOnly && step.momentum().z() <= 0.0) continue;
      if (_firstOnly) {
        auto key = std::make_pair(static_cast<unsigned long>(step.trackId().asUint()), vd);
        auto it = first.find(key);
        if (it == first.end() || step.time() < it->second->time())
          first[key] = &step;
      } else
        kept.push_back(&step);
    }
    for (auto const& kv : first) kept.push_back(kv.second);

    for (const StepPointMC* step : kept) {
      _vdId = static_cast<int>(step->virtualDetectorId());
      _pdgId = step->simParticle()->pdgId();
      _trackId = static_cast<Long64_t>(step->trackId().asUint());
      _x = step->position().x();
      _y = step->position().y();
      _z = step->position().z();
      _px = step->momentum().x();
      _py = step->momentum().y();
      _pz = step->momentum().z();
      _time = step->time();
      double m = 0.0;
      try {
        m = _pdt->particle(_pdgId).mass();
      } catch (const cet::exception&) {
        m = 0.0; // not in the particle table; E is then |p|
      }
      _E = std::sqrt(step->momentum().mag2() + m * m) - m;
      _particles->Fill();
      ++_pdgCounts[_pdgId];
    }
    ++_eventId;
    ++_nEvents;
  };

  void STMParticleDump::endJob() {
    _info->Fill();
    mf::LogInfo log("STMParticleDump");
    log << "read " << _nEvents << " events, wrote " << _particles->GetEntries() << " particles\n";
    for (auto const& kv : _pdgCounts)
      log << "  PDG " << kv.first << ": " << kv.second << "\n";
  };
}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::STMParticleDump)
