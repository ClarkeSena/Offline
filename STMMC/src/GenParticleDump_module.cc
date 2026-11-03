// Writes a GenParticle collection to a flat TTree at full precision, in the same layout as
// STMParticleDump, so countSTMParticleDump.C and the line-window counts work on generator
// output too (e.g. the VDResamplerGenerateMix particles on VD116, tag "VDResamplerGenerator").
// printGenParticle.fcl only prints momenta to 0.1 MeV, too coarse for keV-wide line windows.
//
// Tree "particles", one entry per GenParticle:
//   eventId (L)             sequential input-event index
//   run, subRun, event (I)  art ID of that event
//   virtualdetectorId (I)   fixed value from the fcl (the plane the particles sit on), pdgId (I),
//   trackId (L)             index in the collection
//   x, y, z [mm], px, py, pz [MeV/c], time [ns], E [MeV, kinetic]
// Tree "info", one entry: nEvents (L) = input events read.
// Original author: Sena Clarke

// stdlib includes
#include <cmath>
#include <map>

// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"

// exception handling
#include "cetlib_except/exception.h"

// fhicl includes
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline includes
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"

// ROOT includes
#include "art_root_io/TFileService.h"
#include "TTree.h"

namespace mu2e {
  class GenParticleDump : public art::EDAnalyzer {
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config {
      fhicl::Atom<art::InputTag> GenParticlesTag{Name("GenParticlesTag"), Comment("GenParticle collection, e.g. VDResamplerGenerator")};
      fhicl::Atom<int> virtualDetectorId{Name("virtualDetectorId"), Comment("Value written to the virtualdetectorId branch (the plane the particles are generated on)"), 116};
    };
    using Parameters = art::EDAnalyzer::Table<Config>;
    explicit GenParticleDump(const Parameters& conf);
    void analyze(const art::Event& event) override;
    void endJob() override;

  private:
    art::ProductToken<GenParticleCollection> _genToken;
    int _vdIdConf;
    GlobalConstantsHandle<ParticleDataList> _pdt;

    TTree* _particles = nullptr;
    TTree* _info = nullptr;
    Long64_t _eventId = 0, _trackId = 0, _nEvents = 0;
    int _run = 0, _subRun = 0, _event = 0, _vdId = 0, _pdgId = 0;
    double _x = 0, _y = 0, _z = 0, _px = 0, _py = 0, _pz = 0, _time = 0, _E = 0;
    std::map<int, unsigned long> _pdgCounts;
  };

  GenParticleDump::GenParticleDump(const Parameters& conf) :
    art::EDAnalyzer(conf),
    _genToken(consumes<GenParticleCollection>(conf().GenParticlesTag())),
    _vdIdConf(conf().virtualDetectorId()) {
      art::ServiceHandle<art::TFileService> tfs;
      _particles = tfs->make<TTree>("particles", "Generated particles");
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

  void GenParticleDump::analyze(const art::Event& event) {
    _run = event.run();
    _subRun = event.subRun();
    _event = event.event();
    _vdId = _vdIdConf;
    auto const& gens = event.getProduct(_genToken);
    _trackId = 0;
    for (const GenParticle& gen : gens) {
      _pdgId = static_cast<int>(gen.pdgId());
      _x = gen.position().x();
      _y = gen.position().y();
      _z = gen.position().z();
      _px = gen.momentum().px();
      _py = gen.momentum().py();
      _pz = gen.momentum().pz();
      _time = gen.time();
      double m = 0.0;
      try {
        m = _pdt->particle(_pdgId).mass();
      } catch (const cet::exception&) {
        m = 0.0; // not in the particle table; E is then |p|
      }
      _E = std::sqrt(_px * _px + _py * _py + _pz * _pz + m * m) - m;
      _particles->Fill();
      ++_pdgCounts[_pdgId];
      ++_trackId;
    }
    ++_eventId;
    ++_nEvents;
  };

  void GenParticleDump::endJob() {
    _info->Fill();
    mf::LogInfo log("GenParticleDump");
    log << "GenParticleDump: " << _nEvents << " events, " << _particles->GetEntries() << " particles\n";
    for (auto const& kv : _pdgCounts)
      log << "  pdgId " << kv.first << ": " << kv.second << "\n";
  };
} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::GenParticleDump)
