// Generates photons as a homogeneous Poisson process: the number of photons
// per event is Poisson-distributed and each photon gets an independent uniform
// time in [0, maxTime). Treating each art event as one micropulse window and
// concatenating the digitized waveforms downstream yields a continuous Poisson
// pulse train of rate meanPhotonsPerEvent / maxTime (e.g. 0.0339 / 1695 ns = 20 kHz).
//
// Each photon's ENERGY is drawn from a configurable spectrum:
//   - single fixed line (set E), OR
//   - several discrete lines (lineEnergies) with relative intensities
//     (lineIntensities), plus an optional exponential background continuum
//     (background: { weight, E0, EMin, EMax }), pdf ~ exp(-E/E0) on [EMin,EMax].
// A photon is assigned to a component in proportion to its weight; line photons
// take the line energy, background photons are sampled from the truncated
// exponential by inverse-CDF. The total rate is still meanPhotonsPerEvent; the
// intensities/weight only set the relative mix.
//
// Copied from EventGenerator/src/PhotonGun_module.cc (Claudia Alvarez-Garcia,
// adapted by Pawel Plesniak).
// Original author: Sena Clarke

// stdlib includes
#include <cmath>
#include <limits>
#include <vector>

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// exception handling
#include "cetlib_except/exception.h"

// fhicl includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "fhiclcpp/types/OptionalTable.h"

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline includes
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/SeedService/inc/SeedService.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandFlat.h"

namespace mu2e {
  class PoissonPhotonGun : public art::EDProducer {
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct BackgroundConfig {
      fhicl::Atom<bool> enabled{ Name("enabled"), Comment("Turn the exponential background on/off. Default off."), false};
      fhicl::Atom<double> weight{ Name("weight"), Comment("Relative weight of the exponential background component (same scale as lineIntensities)"), 1.0};
      fhicl::Atom<double> E0{ Name("E0"), Comment("Exponential scale [MeV]; background pdf ~ exp(-E/E0)"), 0.5};
      fhicl::Atom<double> EMin{ Name("EMin"), Comment("Minimum background energy [MeV]"), 0.1};
      fhicl::Atom<double> EMax{ Name("EMax"), Comment("Maximum background energy [MeV]"), 2.0};
    };
    struct Config {
      fhicl::Atom<double> x{ Name("x"), Comment("x position of generated photon")};
      fhicl::Atom<double> y{ Name("y"), Comment("y position of generated photon")};
      fhicl::Atom<double> z{ Name("z"), Comment("z position of generated photon")};
      fhicl::OptionalAtom<double> px{ Name("px"), Comment("x momentum of generated photon")};
      fhicl::OptionalAtom<double> py{ Name("py"), Comment("y momentum of generated photon")};
      fhicl::OptionalAtom<double> E{ Name("E"), Comment("Single-line photon energy [MeV]. Use this OR lineEnergies/background.")};
      fhicl::OptionalSequence<double> lineEnergies{ Name("lineEnergies"), Comment("Discrete photon-line energies [MeV], e.g. [0.347, 0.844, 1.809]")};
      fhicl::OptionalSequence<double> lineIntensities{ Name("lineIntensities"), Comment("Relative intensities for lineEnergies (same length). Default: all 1.0 (equal).")};
      fhicl::OptionalTable<BackgroundConfig> background{ Name("background"), Comment("Optional exponential background continuum")};
      fhicl::Atom<double> meanPhotonsPerEvent{ Name("meanPhotonsPerEvent"), Comment("Poisson mean number of photons per event window. rate = meanPhotonsPerEvent / maxTime")};
      fhicl::Atom<double> maxTime{ Name("maxTime"), Comment("Event window length [ns]; photon times are uniform in [0, maxTime)")};
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit PoissonPhotonGun(const Parameters& conf);
    virtual void produce(art::Event& event);
  private:
    double sampleEnergy();

    double x = 0.0, y = 0.0, z = 0.0;
    double px = 0.0, py = 0.0;
    double meanPhotonsPerEvent = 0.0;
    double maxTime = 0.0;

    // Spectrum: cumulative weights over components. Components are the discrete
    // lines (fixed energy) followed, if present, by the background continuum.
    std::vector<double> _lineE;      // line energies [MeV]
    std::vector<double> _cumWeight;  // cumulative weights, last entry = total
    double _totalWeight = 0.0;
    bool _hasBackground = false;
    double _bgE0 = 0.0, _bgMin = 0.0, _bgMax = 0.0;
    double _bgExpMin = 0.0, _bgExpMax = 0.0; // exp(-EMin/E0), exp(-EMax/E0)

    // Random engine seeded via SeedService, so repeated jobs produce independent
    // Poisson trains (same pattern as the HPGe/LaBr digitizer test modules).
    art::RandomNumberGenerator::base_engine_t& _engine;
    CLHEP::RandPoissonQ _poisson;
    CLHEP::RandFlat _flat;
  };

  PoissonPhotonGun::PoissonPhotonGun(const Parameters& conf):
    art::EDProducer(conf),
    x(conf().x()),
    y(conf().y()),
    z(conf().z()),
    meanPhotonsPerEvent(conf().meanPhotonsPerEvent()),
    maxTime(conf().maxTime()),
    _engine(createEngine(art::ServiceHandle<SeedService>()->getSeed())),
    _poisson(_engine, conf().meanPhotonsPerEvent()),
    _flat(_engine) {
      produces<GenParticleCollection>();
      px = conf().px() ? *conf().px() : 0;
      py = conf().py() ? *conf().py() : 0;
      if (meanPhotonsPerEvent <= 0)
        throw cet::exception("RANGE") << "meanPhotonsPerEvent must be > 0, exiting.";
      if (maxTime <= 0)
        throw cet::exception("RANGE") << "maxTime must be > 0, exiting.";

      // Build the spectrum. Accept either a single E, a set of lines, or both,
      // optionally with an exponential background.
      std::vector<double> intensities;
      std::vector<double> lineEnergies, lineIntensities;
      bool haveLines = conf().lineEnergies(lineEnergies) && !lineEnergies.empty();
      if (haveLines) {
        if (conf().lineIntensities(lineIntensities)) {
          if (lineIntensities.size() != lineEnergies.size())
            throw cet::exception("RANGE") << "lineIntensities must match lineEnergies in length, exiting.";
        } else {
          lineIntensities.assign(lineEnergies.size(), 1.0);
        }
        for (size_t i = 0; i < lineEnergies.size(); ++i) {
          _lineE.push_back(lineEnergies[i]);
          intensities.push_back(lineIntensities[i]);
        }
      }
      double singleE = 0.0;
      if (conf().E(singleE)) {
        _lineE.push_back(singleE);
        intensities.push_back(1.0);
      }

      BackgroundConfig bg;
      if (conf().background(bg) && bg.enabled()) {
        _hasBackground = true;
        _bgE0 = bg.E0();
        _bgMin = bg.EMin();
        _bgMax = bg.EMax();
        if (_bgE0 <= 0.0 || _bgMin < 0.0 || _bgMax <= _bgMin)
          throw cet::exception("RANGE") << "background needs E0>0 and 0<=EMin<EMax, exiting.";
        _bgExpMin = std::exp(-_bgMin / _bgE0);
        _bgExpMax = std::exp(-_bgMax / _bgE0);
        intensities.push_back(bg.weight());
      }

      if (_lineE.empty() && !_hasBackground)
        throw cet::exception("RANGE") << "no spectrum given: set E, lineEnergies, or background, exiting.";

      // Cumulative weights over components (lines first, then background).
      double acc = 0.0;
      for (double w : intensities) {
        if (w < 0.0)
          throw cet::exception("RANGE") << "intensities/weights must be >= 0, exiting.";
        acc += w;
        _cumWeight.push_back(acc);
      }
      _totalWeight = acc;
      if (_totalWeight <= 0.0)
        throw cet::exception("RANGE") << "total spectrum weight must be > 0, exiting.";

      // Guard: transverse momentum must not exceed the smallest possible energy,
      // else pz = sqrt(E^2 - pt^2) is imaginary. Only matters if px/py are set.
      double pt = std::sqrt(px * px + py * py);
      if (pt > 0.0) {
        double eMinPossible = std::numeric_limits<double>::max();
        for (double e : _lineE)
          eMinPossible = std::min(eMinPossible, e);
        if (_hasBackground)
          eMinPossible = std::min(eMinPossible, _bgMin);
        if (pt > eMinPossible)
          throw cet::exception("RANGE") << "transverse momentum " << pt
            << " exceeds the smallest photon energy " << eMinPossible << ", exiting.";
      }
    };

  double PoissonPhotonGun::sampleEnergy() {
    // Pick a component in proportion to its weight.
    double u = _flat.fire() * _totalWeight;
    size_t i = 0;
    for (; i < _cumWeight.size(); ++i)
      if (u <= _cumWeight[i])
        break;
    if (i < _lineE.size())
      return _lineE[i];
    // Background component: inverse-CDF sample of exp(-E/E0) truncated to [EMin,EMax].
    double U = _flat.fire();
    return -_bgE0 * std::log(_bgExpMin - U * (_bgExpMin - _bgExpMax));
  }

  void PoissonPhotonGun::produce(art::Event& event) {
    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);
    const CLHEP::Hep3Vector pos(x, y, z);
    int nPhotons = _poisson.fire();
    for (int i = 0; i < nPhotons; i++) {
      double E = sampleEnergy();
      double pz = std::sqrt(std::max(0.0, E * E - px * px - py * py));
      const CLHEP::Hep3Vector p(px, py, pz);
      CLHEP::HepLorentzVector mom(p, E);
      double t = _flat.fire(maxTime);
      output->push_back(GenParticle(PDGCode::gamma, GenId::particleGun, pos, mom, t));
    }
    event.put(std::move(output));
  };
}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::PoissonPhotonGun)
