// Monoenergetic photon line emitted from a muon stop, e.g. the muonic Al 2p-1s X-ray (346.8 keV)
// or the 844 keV gamma from 27Mg beta decay following capture. Generalisation of
// MuCap1809keVGammaGenerator with the energy, intensity and emission delay set in fhicl, so a new
// line needs no new code or global constants.
//
// The intensity is per call of generate(), so it must match the Pileup list the tool is placed in:
// per stop for stopProducts, per capture for captureProducts, per decay for decayProducts.
//
// Emission delay (Kinematic::timeDelay, relative to the time the Pileup module hands over):
//   lifetime = 0                 : prompt
//   lifetime > 0, foldPeriod = 0 : Exp(lifetime)
//   lifetime > 0, foldPeriod > 0 : Exp(lifetime) folded into [0, foldPeriod). For a long-lived
//                                  emitter (lifetime >> foldPeriod) this is uniform in time, i.e. a
//                                  constant rate uncorrelated with the bunch structure, without
//                                  pushing photon times out to seconds.

#include "art/Utilities/ToolMacros.h"
#include "cetlib_except/exception.h"

#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Random/RandFlat.h"

#include "Offline/EventGenerator/inc/ParticleGeneratorTool.hh"

#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/Mu2eUtilities/inc/RandomUnitSphere.hh"

#include <cmath>

namespace mu2e {
  class MuStopGammaLineGenerator : public ParticleGeneratorTool {
  public:
    struct PhysConfig {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<double> energy     {Name("energy")    , Comment("Photon energy (MeV)")};
      fhicl::Atom<double> intensity  {Name("intensity") , Comment("Emission probability per generate() call, in [0, 1]")};
      fhicl::Atom<double> lifetime   {Name("lifetime")  , Comment("Mean emission delay (ns), 0 for prompt"), 0.};
      fhicl::Atom<double> foldPeriod {Name("foldPeriod"), Comment("If > 0, fold the delay into [0, foldPeriod) (ns)"), 0.};
      fhicl::Atom<double> czMin      {Name("czmin")     , Comment("Restrict cos(theta_z) minimum"), -1.};
      fhicl::Atom<double> czMax      {Name("czmax")     , Comment("Restrict cos(theta_z) maximum"),  1.};
    };
    typedef art::ToolConfigTable<PhysConfig> Parameters;

    explicit MuStopGammaLineGenerator(Parameters const& conf) :
      _energy(conf().energy()),
      _intensity(conf().intensity()),
      _lifetime(conf().lifetime()),
      _foldPeriod(conf().foldPeriod()),
      _czMin(conf().czMin()),
      _czMax(conf().czMax()),
      _pdgId(PDGCode::gamma)
    {
      if(_energy <= 0.) throw cet::exception("BADCONFIG") << "Photon energy must be positive: " << _energy;
      if(_intensity < 0. || _intensity > 1.) throw cet::exception("BADCONFIG") << "Intensity outside [0, 1]: " << _intensity;
      if(_lifetime < 0. || _foldPeriod < 0.) throw cet::exception("BADCONFIG") << "lifetime and foldPeriod must be >= 0";
      if(_czMin < -1. || _czMax > 1. || _czMin > _czMax) throw cet::exception("BADCONFIG") << "Cos(theta_z) range unphysical: " << _czMin << " - " << _czMax;
    }

    std::vector<ParticleGeneratorTool::Kinematic> generate() override;
    void generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) override;
    std::unique_ptr<SpectrumConfig> spectrumConfig() override;

    void finishInitialization(art::RandomNumberGenerator::base_engine_t& eng, const std::string& /*material*/, const bool isPrimary) override {
      _isPrimary = isPrimary;
      _randomUnitSphere = std::make_unique<RandomUnitSphere>(eng, _czMin, _czMax);
      _randFlat = std::make_unique<CLHEP::RandFlat>(eng);
      _randExp = std::make_unique<CLHEP::RandExponential>(eng);
    }

  private:
    double _energy;
    double _intensity;
    double _lifetime;
    double _foldPeriod;
    double _czMin;
    double _czMax;
    PDGCode::type _pdgId;

    std::unique_ptr<RandomUnitSphere>       _randomUnitSphere;
    std::unique_ptr<CLHEP::RandFlat>        _randFlat;
    std::unique_ptr<CLHEP::RandExponential> _randExp;

    double emissionDelay();
  };

  double MuStopGammaLineGenerator::emissionDelay() {
    if(_lifetime <= 0.) return 0.;
    const double delay = _randExp->fire(_lifetime);
    return (_foldPeriod > 0.) ? std::fmod(delay, _foldPeriod) : delay;
  }

  std::vector<ParticleGeneratorTool::Kinematic> MuStopGammaLineGenerator::generate() {
    std::vector<ParticleGeneratorTool::Kinematic>  res;

    const double intensity = _intensity * (_czMax - _czMin)/2.; // account for potential cz selection in the intensity
    if (_isPrimary || _randFlat->fire() < intensity) {
      CLHEP::Hep3Vector p3 = _randomUnitSphere->fire(_energy); // massless: |p| = E
      CLHEP::HepLorentzVector fourmom(p3, _energy);
      ParticleGeneratorTool::Kinematic k{_pdgId, ProcessCode::mu2eMuonCaptureAtRest, fourmom};
      k.timeDelay = emissionDelay();
      res.emplace_back(k);
    }

    return res;
  }

  void MuStopGammaLineGenerator::generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) {
    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);
    const auto daughters = generate();
    for(const auto& d: daughters) {
      out->emplace_back(d.pdgId,
                        GenId::MuCapGammaRayGenTool,
                        pos,
                        d.fourmom,
                        stop.t + d.timeDelay);
    }
  }

  std::unique_ptr<SpectrumConfig> MuStopGammaLineGenerator::spectrumConfig() {
    auto config = std::make_unique<SpectrumConfig>();
    config->add_var(SpectrumConfig::RestrictedVar("energy", 1., _energy, _energy));
    config->add_var(SpectrumConfig::RestrictedVar("cosz", (_czMax - _czMin)/2., _czMin, _czMax));
    return config;
  }

}
DEFINE_ART_CLASS_TOOL(mu2e::MuStopGammaLineGenerator)
