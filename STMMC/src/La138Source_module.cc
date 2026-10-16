// Intrinsic LaBr3 radioactivity source, La-138 only (Method C): for each La-138 decay, drawn
// at the rate implied by the crystal's lanthanum content, emit the decay PRODUCTS directly as
// GenParticles at a uniform vertex inside the crystal. Geant4 then transports them, so the
// photopeak/Compton/escape structure, the EC coincidence sum (1436 keV gamma + Ba atomic
// relaxation = 1473 keV) and the beta- 789 keV-plus-continuum emerge from transport. No G4
// radioactive decay is used (set turnOnRadioactiveDecay: false): analogue G4RadioactiveDecay
// never decays La-138 in-window because its 1e11 yr half-life gives a ~1e-25 probability.
//
// Decay kinematics live in Offline/STMMC/inc/LaBrIntrinsicDecays.hh (shared with
// LaBrIntrinsicSource, which adds the Ac-227 chain). Relative to the first version:
//   - EC atomic relaxation is complete (K X-ray + remaining cascade energy, all-Auger K,
//     L and M capture, with DDEP capture ratios), so K capture sums to 1473 keV and L
//     capture to 1442 keV, not 1468 or 1436.
//   - The beta spectrum uses the 2nd-forbidden-unique shape factor and a Fermi function.
//   - Ce (and Sr) substitution on La sites is included via laPerMolecule (default 0.95).
//   - The activity can instead be fixed to a literature specific activity [Bq/cm^3].
//
// Activity. By default COMPUTED from the geometry and composition:
//   A = ln2 / T1/2 * (density * V / molarMass) * N_A * laPerMolecule * abundance(La-138)
// With the defaults this gives ~1.454 Bq/cm^3, i.e. ~505 Bq for a 3"x3" crystal. Literature:
//   1.425(59) Bq/cm^3  Saint-Gobain B380 3"x3" (Cheng et al., NST 31 (2020) 99, arXiv:2007.07552)
//   1.488 Bq/cm^3      average over several LaBr3(Ce) crystals quoted in the same paper
// Setting la138SpecificActivity overrides the computed value (the computed one is still
// printed for comparison). Either way the rate scales with the crystal VOLUME given here.
//
// Rate per event: mean decays = A[Bq] * activityScale * windowTime[ns] * 1e-9. Products are
// tagged GenId::STMLa138 so the truth modules can separate them from other sources.
// Original author: Sena Clarke

// stdlib includes
#include <cmath>
#include <string>

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// exception handling
#include "cetlib_except/exception.h"

// fhicl includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"

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
  class La138Source : public art::EDProducer {
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
      fhicl::Atom<double> windowTime{ Name("windowTime"), Comment("Event time window [ns] over which to integrate the rate; match micropulseTime.")};

      // Activity: literature value, or computed from the composition below.
      fhicl::OptionalAtom<double> la138SpecificActivity{ Name("la138SpecificActivity"), Comment("If set, use this La-138 specific activity [Bq/cm^3] instead of the computed one. Literature: 1.425 (B380 3x3, Cheng 2020), 1.488 (multi-crystal average).")};
      fhicl::OptionalAtom<double> density{ Name("density"), Comment("Crystal density [g/cm^3]. Default 5.08 (LaBr3:Ce).")};
      fhicl::OptionalAtom<double> molarMass{ Name("molarMass"), Comment("Molar mass [g/mol]. Default 378.62 (LaBr3).")};
      fhicl::OptionalAtom<double> laPerMolecule{ Name("laPerMolecule"), Comment("La atoms per formula unit; 1 - dopant fraction. Default 0.95 (5% Ce).")};
      fhicl::OptionalAtom<double> la138Abundance{ Name("la138Abundance"), Comment("La-138 isotopic abundance. Default 8.881e-4 (IUPAC).")};
      fhicl::OptionalAtom<double> halfLife_yr{ Name("halfLife_yr"), Comment("La-138 half-life [yr]. Default 1.03e11 (ENSDF/NUBASE2020).")};

      // Decay scheme and model options.
      fhicl::OptionalAtom<double> ecFraction{ Name("ecFraction"), Comment("EC branch fraction. Default 0.655 (ENSDF); DDEP 0.652.")};
      fhicl::OptionalAtom<std::string> betaShape{ Name("betaShape"), Comment("'secondForbiddenUnique' (default, correct for La-138) or 'allowed' (previous behaviour).")};

      fhicl::OptionalAtom<double> activityScale{ Name("activityScale"), Comment("Multiplier on the activity for statistics studies. Default 1 (physical).")};
      fhicl::OptionalAtom<int> verbosityLevel{ Name("verbosityLevel"), Comment("Verbosity. Default 0.")};
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit La138Source(const Parameters& conf);
    void produce(art::Event& event) override;
  private:
    static LaBrIntrinsic::BetaShape parseShape(const std::string& s);

    CLHEP::Hep3Vector centre;
    double radius = 0, halfLength = 0;
    double windowTime = 0;
    double meanDecaysPerEvent = 0;
    int verbosityLevel = 0;

    art::RandomNumberGenerator::base_engine_t& _engine;
    CLHEP::RandFlat _flat;
    CLHEP::RandPoissonQ _poisson;
    LaBrIntrinsic::La138Decay _decay;
  };

  LaBrIntrinsic::BetaShape La138Source::parseShape(const std::string& s) {
    if (s == "secondForbiddenUnique") return LaBrIntrinsic::BetaShape::secondForbiddenUnique;
    if (s == "allowed")               return LaBrIntrinsic::BetaShape::allowed;
    throw cet::exception("RANGE") << "betaShape must be 'secondForbiddenUnique' or 'allowed'\n";
  }

  La138Source::La138Source(const Parameters& conf):
    art::EDProducer(conf),
    centre(conf().x(), conf().y(), conf().z()),
    radius(conf().radius()), halfLength(conf().halfLength()),
    windowTime(conf().windowTime()),
    _engine(createEngine(art::ServiceHandle<SeedService>()->getSeed())),
    _flat(_engine),
    _poisson(_engine, 1.0),
    _decay([&] {
      LaBrIntrinsic::La138Data d;
      if (conf().ecFraction()) d.ecFraction = *(conf().ecFraction());
      if (d.ecFraction < 0.0 || d.ecFraction > 1.0)
        throw cet::exception("RANGE") << "ecFraction must be in [0,1]\n";
      return d;
    }(), parseShape(conf().betaShape() ? *(conf().betaShape()) : "secondForbiddenUnique")) {
      produces<GenParticleCollection>();

      const double density       = conf().density()        ? *(conf().density())        : 5.08;
      const double molarMass     = conf().molarMass()      ? *(conf().molarMass())      : 378.62;
      const double laPerMolecule = conf().laPerMolecule()  ? *(conf().laPerMolecule())  : 0.95;
      const double la138Abund    = conf().la138Abundance() ? *(conf().la138Abundance()) : 8.881e-4;
      const double halfLife_yr   = conf().halfLife_yr()    ? *(conf().halfLife_yr())    : 1.03e11;
      const double activityScale = conf().activityScale()  ? *(conf().activityScale())  : 1.0;
      verbosityLevel             = conf().verbosityLevel() ? *(conf().verbosityLevel()) : 0;

      if (radius <= 0 || halfLength <= 0 || windowTime <= 0)
        throw cet::exception("RANGE") << "radius, halfLength and windowTime must be > 0\n";
      if (activityScale < 0)
        throw cet::exception("RANGE") << "activityScale must be >= 0\n";

      const double V_cm3 = LaBrIntrinsic::cylinderVolume_cm3(radius, halfLength);
      const double computed = LaBrIntrinsic::la138SpecificActivity(
        density, molarMass, laPerMolecule, la138Abund, halfLife_yr);
      const bool useLiterature = static_cast<bool>(conf().la138SpecificActivity());
      const double specific = useLiterature ? *(conf().la138SpecificActivity()) : computed;
      if (specific < 0)
        throw cet::exception("RANGE") << "la138SpecificActivity must be >= 0\n";

      const double activity_Bq = specific * V_cm3 * activityScale;
      meanDecaysPerEvent = activity_Bq * windowTime * 1e-9;   // window ns -> s

      if (verbosityLevel > 0) {
        mf::LogInfo("La138Source")
          << "crystal V = " << V_cm3 << " cm^3, mass = " << density * V_cm3 << " g\n"
          << "  computed specific activity = " << computed << " Bq/cm^3"
          << " (laPerMolecule " << laPerMolecule << ")\n"
          << "  literature: 1.425(59) Bq/cm^3 (B380 3x3), 1.488 Bq/cm^3 (multi-crystal avg)\n"
          << "  USING " << specific << " Bq/cm^3 (" << (useLiterature ? "fcl" : "computed") << ")"
          << " -> activity = " << activity_Bq << " Bq" << (activityScale != 1.0 ? " (scaled)" : "") << "\n"
          << "  EC fraction = " << _decay.data().ecFraction << "\n"
          << "  mean decays per " << windowTime << " ns event = " << meanDecaysPerEvent;
      }
    }

  void La138Source::produce(art::Event& event) {
    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);
    const int nDecays = _poisson.fire(meanDecaysPerEvent);
    for (int i = 0; i < nDecays; ++i) {
      LaBrIntrinsic::Emitter e{_flat, *output, GenId::STMLa138,
                               LaBrIntrinsic::sampleCylinder(_flat, centre, radius, halfLength),
                               _flat.fire() * windowTime};   // uniform in window: constant activity
      _decay.emit(e);
    }
    event.put(std::move(output));
  }
} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::La138Source)
