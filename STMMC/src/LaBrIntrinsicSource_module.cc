// Full intrinsic radioactivity of a LaBr3 crystal: La-138 PLUS the Ac-227 decay chain.
//
// Same method as La138Source (decay products emitted directly at uniform vertices inside the
// crystal, no G4 radioactive decay), with the same La-138 physics from
// Offline/STMMC/inc/LaBrIntrinsicDecays.hh. On top of that, each of the eight Ac-227 chain
// members (Ac-227, Th-227, Ra-223, Rn-219, Po-215, Pb-211, Bi-211, Tl-207) is an independent
// Poisson process at its own activity.
//
// Why the Ac-227 chain matters for Mu2e: its 5.7-7.4 MeV alphas sit near the 1809 keV signal
// line once quenched. In the literature (alpha/electron ~0.3) they land at ~1.5-2.4 MeV,
// directly under it; in the Mu2e crystal (bench-measured ~0.46) they land at 2.5-3.7 MeV,
// above it. Contamination varies by crystal (up to ~100x), so it must be fitted to bench
// data. Alphas deposit their FULL energy in Geant4; set alphaQuenchingFactor on the digitizer
// and truth modules (see LaBrTruthHelpers.hh) or the alphas appear at 5-7.4 MeV.
//
// Ac-227 activity (Bq/cm^3, per chain member in secular equilibrium):
//   Mu2e LaBr3:Ce+Sr crystal, bench Eu-152 run 125382 (2026-10-19, labrSelfActivity.py):
//     17.2 +- 0.4 alphas/s in the 2.4-3.8 MeV band -> default ac227ChainActivity 0.0099,
//     total alpha 0.0495 Bq/cm^3 (5 alpha emitters), ~1/30 of La-138. See prolog.fcl.
//   Literature, Saint-Gobain B380 3"x3" LaBr3:Ce (Cheng et al., NST 31 (2020) 99,
//   arXiv:2007.07552): Th-227 0.0158(22), Ra-223 0.0127(14), Rn-219 0.0125(17),
//   Bi-211 0.0136(15) -> 0.0136 (the previous default), Th227+Ra223+Rn219+Bi211 0.0545.
//   Multi-crystal average quoted in the same paper: 0.044 Bq/cm^3 total alpha.
// In equilibrium each member has the chain activity times its branching fraction (Th-227
// 0.9862, Tl-207 0.9972, the rest 1). memberSpecificActivities overrides this per member,
// e.g. to fit a crystal out of equilibrium.
//
// Products are tagged GenId::STMLa138 or GenId::STMAc227 so the truth modules report the two
// separately. Rate per event for each species: A[Bq] * activityScale * windowTime[ns] * 1e-9.
// Original author: Sena Clarke

// stdlib includes
#include <cmath>
#include <sstream>
#include <string>
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

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline includes
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/SeedService/inc/SeedService.hh"
#include "Offline/STMMC/inc/LaBrIntrinsicDecays.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandFlat.h"

namespace mu2e {
  class LaBrIntrinsicSource : public art::EDProducer {
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config {
      // Crystal placement and size (Mu2e global coords, mm).
      fhicl::Atom<double> x{ Name("x"), Comment("Crystal centre x [mm]")};
      fhicl::Atom<double> y{ Name("y"), Comment("Crystal centre y [mm]")};
      fhicl::Atom<double> z{ Name("z"), Comment("Crystal centre z [mm]")};
      fhicl::Atom<double> radius{ Name("radius"), Comment("Crystal radius [mm]")};
      fhicl::Atom<double> halfLength{ Name("halfLength"), Comment("Crystal half-length along z [mm]")};
      fhicl::Atom<double> windowTime{ Name("windowTime"), Comment("Event time window [ns]; match micropulseTime.")};

      // ---- La-138 (identical meaning to La138Source) ----
      fhicl::Atom<bool> includeLa138{ Name("includeLa138"), Comment("Generate La-138 decays. Default true; false gives the Ac-227 chain alone."), true};
      fhicl::OptionalAtom<double> la138SpecificActivity{ Name("la138SpecificActivity"), Comment("If set, La-138 specific activity [Bq/cm^3] instead of the computed one.")};
      fhicl::OptionalAtom<double> density{ Name("density"), Comment("Crystal density [g/cm^3]. Default 5.08.")};
      fhicl::OptionalAtom<double> molarMass{ Name("molarMass"), Comment("Molar mass [g/mol]. Default 378.62.")};
      fhicl::OptionalAtom<double> laPerMolecule{ Name("laPerMolecule"), Comment("La atoms per formula unit. Default 0.95 (5% Ce).")};
      fhicl::OptionalAtom<double> la138Abundance{ Name("la138Abundance"), Comment("La-138 isotopic abundance. Default 8.881e-4.")};
      fhicl::OptionalAtom<double> halfLife_yr{ Name("halfLife_yr"), Comment("La-138 half-life [yr]. Default 1.03e11 (ENSDF).")};
      fhicl::OptionalAtom<double> ecFraction{ Name("ecFraction"), Comment("La-138 EC branch fraction. Default 0.655 (ENSDF).")};

      // ---- Ac-227 chain ----
      fhicl::Atom<bool> includeAc227{ Name("includeAc227"), Comment("Generate the Ac-227 chain. Default true."), true};
      fhicl::Atom<double> ac227ChainActivity{ Name("ac227ChainActivity"), Comment("Ac-227 specific activity [Bq/cm^3]; daughters follow in secular equilibrium. Default 0.0099 (Mu2e crystal, bench run 125382; literature B380 3x3, Cheng 2020: 0.0136)."), 0.0099};
      fhicl::OptionalSequence<double> memberSpecificActivities{ Name("memberSpecificActivities"), Comment("Override: 8 specific activities [Bq/cm^3] in order Ac227, Th227, Ra223, Rn219, Po215, Pb211, Bi211, Tl207.")};

      fhicl::OptionalAtom<double> activityScale{ Name("activityScale"), Comment("Multiplier on ALL activities, for statistics. Default 1.")};
      fhicl::OptionalAtom<int> verbosityLevel{ Name("verbosityLevel"), Comment("Verbosity. Default 0.")};
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit LaBrIntrinsicSource(const Parameters& conf);
    void produce(art::Event& event) override;
    void endJob() override;
  private:
    CLHEP::Hep3Vector centre;
    double radius = 0, halfLength = 0;
    double windowTime = 0;
    int verbosityLevel = 0;

    art::RandomNumberGenerator::base_engine_t& _engine;
    CLHEP::RandFlat _flat;
    CLHEP::RandPoissonQ _poisson;

    LaBrIntrinsic::La138Decay _la138;
    double _la138Mean = 0.0;                        // mean La-138 decays per event
    std::vector<LaBrIntrinsic::Nuclide> _chain;
    std::vector<double> _chainMean;                 // mean decays per event, per member

    // Diagnostics
    uint64_t n_events = 0, n_la138 = 0;
    std::vector<uint64_t> n_chain;
  };

  LaBrIntrinsicSource::LaBrIntrinsicSource(const Parameters& conf):
    art::EDProducer(conf),
    centre(conf().x(), conf().y(), conf().z()),
    radius(conf().radius()), halfLength(conf().halfLength()),
    windowTime(conf().windowTime()),
    _engine(createEngine(art::ServiceHandle<SeedService>()->getSeed())),
    _flat(_engine),
    _poisson(_engine, 1.0),
    _la138([&] {
      LaBrIntrinsic::La138Data d;
      if (conf().ecFraction()) d.ecFraction = *(conf().ecFraction());
      if (d.ecFraction < 0.0 || d.ecFraction > 1.0)
        throw cet::exception("RANGE") << "ecFraction must be in [0,1]\n";
      return d;
    }()),
    _chain(LaBrIntrinsic::makeAc227Chain()) {
      produces<GenParticleCollection>();

      const double activityScale = conf().activityScale()  ? *(conf().activityScale())  : 1.0;
      verbosityLevel             = conf().verbosityLevel() ? *(conf().verbosityLevel()) : 0;
      if (radius <= 0 || halfLength <= 0 || windowTime <= 0)
        throw cet::exception("RANGE") << "radius, halfLength and windowTime must be > 0\n";
      if (activityScale < 0)
        throw cet::exception("RANGE") << "activityScale must be >= 0\n";

      const double V_cm3 = LaBrIntrinsic::cylinderVolume_cm3(radius, halfLength);
      const double toMeanPerEvent = V_cm3 * activityScale * windowTime * 1e-9;

      // ---- La-138 ----
      const double density       = conf().density()        ? *(conf().density())        : 5.08;
      const double molarMass     = conf().molarMass()      ? *(conf().molarMass())      : 378.62;
      const double laPerMolecule = conf().laPerMolecule()  ? *(conf().laPerMolecule())  : 0.95;
      const double la138Abund    = conf().la138Abundance() ? *(conf().la138Abundance()) : 8.881e-4;
      const double halfLife_yr   = conf().halfLife_yr()    ? *(conf().halfLife_yr())    : 1.03e11;
      const double laComputed = LaBrIntrinsic::la138SpecificActivity(
        density, molarMass, laPerMolecule, la138Abund, halfLife_yr);
      const bool laFromFcl = static_cast<bool>(conf().la138SpecificActivity());
      const double laSpecific = !conf().includeLa138() ? 0.0
                              : laFromFcl ? *(conf().la138SpecificActivity()) : laComputed;
      if (laSpecific < 0)
        throw cet::exception("RANGE") << "la138SpecificActivity must be >= 0\n";
      _la138Mean = laSpecific * toMeanPerEvent;

      // ---- Ac-227 chain ----
      std::vector<double> specific(_chain.size(), 0.0);
      std::vector<double> overrides;
      const bool haveOverrides = conf().memberSpecificActivities(overrides);
      if (haveOverrides && overrides.size() != _chain.size())
        throw cet::exception("RANGE") << "memberSpecificActivities needs exactly "
                                      << _chain.size() << " values\n";
      if (conf().ac227ChainActivity() < 0)
        throw cet::exception("RANGE") << "ac227ChainActivity must be >= 0\n";
      for (size_t i = 0; i < _chain.size(); ++i) {
        if (!conf().includeAc227()) break;
        specific[i] = haveOverrides ? overrides[i]
                                    : conf().ac227ChainActivity() * _chain[i].chainFraction();
        if (specific[i] < 0)
          throw cet::exception("RANGE") << "memberSpecificActivities must be >= 0\n";
      }
      for (double s : specific) _chainMean.push_back(s * toMeanPerEvent);
      n_chain.assign(_chain.size(), 0);

      if (verbosityLevel > 0) {
        std::ostringstream os;
        double chainTotal = 0.0, alphaTotal = 0.0;
        for (size_t i = 0; i < _chain.size(); ++i) {
          chainTotal += specific[i];
          // Th-227, Ra-223, Rn-219, Bi-211: the alpha emitters Cheng et al. quantified via
          // their gamma lines, so this sum is directly comparable to their 0.0545 Bq/cm^3
          // (and to 4 x 0.0099 = 0.0396 Bq/cm^3 measured on the Mu2e crystal).
          if (i == 1 || i == 2 || i == 3 || i == 6) alphaTotal += specific[i];
          os << "\n    " << _chain[i].name() << " : " << specific[i] << " Bq/cm^3 -> "
             << specific[i] * V_cm3 << " Bq";
        }
        mf::LogInfo("LaBrIntrinsicSource")
          << "crystal V = " << V_cm3 << " cm^3\n"
          << "  La-138: computed " << laComputed << " Bq/cm^3, USING " << laSpecific
          << " Bq/cm^3 (" << (!conf().includeLa138() ? "disabled" : laFromFcl ? "fcl" : "computed")
          << ") -> " << laSpecific * V_cm3 << " Bq\n"
          << "  Ac-227 chain (" << (!conf().includeAc227() ? "disabled"
                                    : haveOverrides ? "per-member fcl" : "secular equilibrium") << "):"
          << os.str() << "\n"
          << "    total " << chainTotal << " Bq/cm^3; Th227+Ra223+Rn219+Bi211 = " << alphaTotal
          << " Bq/cm^3 (Mu2e bench: 0.0396; literature B380: 0.0545)\n"
          << "  La-138 / (Th227+Ra223+Rn219+Bi211) = " << (alphaTotal > 0 ? laSpecific / alphaTotal : 0.0)
          << "  (Mu2e bench: 37; literature B380: ~25)\n"
          << "  activityScale = " << activityScale
          << ", mean La-138 decays per event = " << _la138Mean;
      }
    }

  void LaBrIntrinsicSource::produce(art::Event& event) {
    ++n_events;
    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);

    auto newEmitter = [&](GenId id) {
      return LaBrIntrinsic::Emitter{_flat, *output, id,
        LaBrIntrinsic::sampleCylinder(_flat, centre, radius, halfLength),
        _flat.fire() * windowTime};
    };

    if (_la138Mean > 0) {
      const int n = _poisson.fire(_la138Mean);
      for (int i = 0; i < n; ++i) {
        auto e = newEmitter(GenId::STMLa138);
        _la138.emit(e);
      }
      n_la138 += n;
    }
    for (size_t k = 0; k < _chain.size(); ++k) {
      if (_chainMean[k] <= 0) continue;
      const int n = _poisson.fire(_chainMean[k]);
      for (int i = 0; i < n; ++i) {
        auto e = newEmitter(GenId::STMAc227);
        _chain[k].emit(e);
      }
      n_chain[k] += n;
    }
    event.put(std::move(output));
  }

  void LaBrIntrinsicSource::endJob() {
    std::ostringstream os;
    for (size_t k = 0; k < _chain.size(); ++k)
      os << " " << _chain[k].name() << "=" << n_chain[k];
    mf::LogInfo("LaBrIntrinsicSource")
      << "events: " << n_events << ", La-138 decays: " << n_la138
      << ", Ac-227 chain decays:" << os.str();
  }
} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::LaBrIntrinsicSource)
