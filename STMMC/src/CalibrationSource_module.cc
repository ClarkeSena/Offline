// Sealed radioactive calibration source placed at a distance from the LaBr detector.
//
// Emits the gamma rays of a chosen isotope from a point vertex at the ABSOLUTE rate implied
// by a stated source activity, using the same rate convention as La138Source_module.cc:
//
//   meanDecaysPerEvent = activity[Bq] * activityScale * windowTime[ns] * 1e-9
//
// Because both generators express their rate the same way, running this module together with
// La138Source in one job puts a source of KNOWN activity and the crystal's intrinsic La-138
// activity on the same absolute footing. The ratio of their count rates is then a prediction
// that can be checked against a bench measurement, which is what lets the La-138 rate be
// validated instead of scaled by hand.
//
// Directional biasing. A point source at a distance emits isotropically, so only the small
// solid angle subtended by the crystal is useful; at 200 mm that is under 1%, and simulating
// the other 99% is wasted Geant4 time. With emissionMode "cone", each gamma is independently
// kept with probability f = (1-cos(thetaMax))/2 and, if kept, given a direction uniform inside
// a cone of half-angle thetaMax about the source-to-crystal axis. That is mathematically
// identical to isotropic emission for every particle that could reach the crystal, with unit
// weights throughout, and it costs a factor ~1/f less CPU.
//
// The cone is sized from a sphere bounding the whole crystal, so it is valid for any source
// position, not just on-axis. What biasing DOES discard is a gamma emitted away from the
// detector that scatters off surrounding material back into it. On a bench that is a small
// correction; inside the full Mu2e geometry it is not necessarily, so cross-check any final
// number with emissionMode "isotropic".
//
// Placement and rate, for matching a bench measurement:
//   - x, y, z give the source position. With coordinates: "crystal" they are an offset from
//     the crystal centre (crystal axis = z, front face at z = -crystalHalfLength), so a bench
//     layout can be entered directly, e.g. (0, 0, -138.1) is on-axis 100 mm from the face,
//     (60, 0, 0) is 60 mm from the axis beside the crystal. With "global" (default) they are
//     Mu2e coordinates. distanceFromCrystalFace remains as an on-axis shortcut.
//   - activity [Bq] or activity_uCi. If the certificate activity is at a reference date, give
//     elapsedDays since then and the activity is decay-corrected with the isotope's half-life
//     (halfLife_days is required for 'custom').
// The vertex must lie in air/vacuum in the geometry being used: inside the Mu2e STM geometry,
// check that the chosen point is not inside shielding (a point in lead is silently absorbed).
//
// Decay data are ENSDF values (IAEA LiveChart of Nuclides), checked 2026-09-24.
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
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/SeedService/inc/SeedService.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandFlat.h"

namespace {
  // ---- Decay data: ENSDF via the IAEA LiveChart of Nuclides (checked 2026-09-24). ----
  // Energies in MeV, probabilities per DECAY.

  // 137Cs -> 137mBa -> 137Ba. The 662 keV gamma is emitted in 85.1% of decays; the
  // remainder proceed by internal conversion, which leaves a K vacancy in the Ba
  // daughter and yields K X-rays in ~7.1% of all decays.
  constexpr double kCs137_Gamma      = 0.661657;
  constexpr double kCs137_GammaProb  = 0.851;
  constexpr double kCs137_KaE        = 0.032060;  // Ba K-alpha (intensity-weighted)
  constexpr double kCs137_KaProb     = 0.0566;   // Ka2 1.99% + Ka1 3.67%
  constexpr double kCs137_KbE        = 0.036400;  // Ba K-beta
  constexpr double kCs137_KbProb     = 0.0135;   // total K-beta

  // 60Co -> 60Ni*. Two gammas in a TRUE CASCADE: the intermediate state lives ~0.7 ps,
  // so both are emitted from the same vertex at the same time. That is what produces the
  // 2505 keV coincidence sum peak when both are absorbed.
  constexpr double kCo60_Gamma1      = 1.173228;
  constexpr double kCo60_Gamma1Prob  = 0.9985;
  constexpr double kCo60_Gamma2      = 1.332492;
  constexpr double kCo60_Gamma2Prob  = 0.999826;

  // 22Na -> 22Ne*. 1275 keV gamma, plus positron annihilation quanta from the 90.3%
  // beta+ branch. The positron is assumed to annihilate at rest inside the source
  // encapsulation, which is the normal situation for a sealed calibration source.
  constexpr double kNa22_Gamma       = 1.274537;
  constexpr double kNa22_GammaProb   = 0.9994;
  constexpr double kNa22_BetaPlusBr  = 0.899;
  constexpr double kNa22_Annih       = 0.510999;

  // 54Mn -> 54Cr*. Single line. The Cr K X-rays (~5.4 keV) are far below any STM
  // threshold and are not emitted.
  constexpr double kMn54_Gamma       = 0.834848;
  constexpr double kMn54_GammaProb   = 0.99976;
}

namespace mu2e {
  class CalibrationSource : public art::EDProducer {
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config {
      fhicl::Atom<std::string> isotope{ Name("isotope"), Comment("Cs137, Co60, Na22, Mn54, or custom (then set lineEnergies/lineIntensities)")};
      fhicl::OptionalAtom<double> activity{ Name("activity"), Comment("Source activity [Bq]. Give this OR activity_uCi.")};
      fhicl::OptionalAtom<double> activity_uCi{ Name("activity_uCi"), Comment("Source activity [uCi] (1 uCi = 37000 Bq). Give this OR activity.")};
      fhicl::OptionalAtom<double> elapsedDays{ Name("elapsedDays"), Comment("Days from the activity reference date to the measurement; the activity is decay-corrected. Default 0.")};
      fhicl::OptionalAtom<double> halfLife_days{ Name("halfLife_days"), Comment("Half-life [days] for decay correction. Built in for Cs137/Co60/Na22/Mn54; required for 'custom' if elapsedDays > 0.")};
      fhicl::Atom<double> windowTime{ Name("windowTime"), Comment("Event time window [ns] over which to integrate the rate; match micropulseTime.")};
      fhicl::OptionalAtom<double> activityScale{ Name("activityScale"), Comment("Multiplier on the activity for statistics studies. Set the SAME value on La138Source to preserve the rate ratio. Default 1 (physical).")};

      // Custom spectrum: independent lines, each emitted with its own per-decay probability.
      fhicl::OptionalSequence<double> lineEnergies{ Name("lineEnergies"), Comment("Custom gamma energies [MeV]; used when isotope is 'custom'.")};
      fhicl::OptionalSequence<double> lineIntensities{ Name("lineIntensities"), Comment("Emission probability PER DECAY for each custom line (absolute, not relative). Same length as lineEnergies.")};

      // Crystal, used both for on-axis placement and for sizing the bias cone.
      fhicl::OptionalAtom<double> crystalX{ Name("crystalX"), Comment("LaBr crystal centre x [mm]. Default -3863.36.")};
      fhicl::OptionalAtom<double> crystalY{ Name("crystalY"), Comment("LaBr crystal centre y [mm]. Default 0.")};
      fhicl::OptionalAtom<double> crystalZ{ Name("crystalZ"), Comment("LaBr crystal centre z [mm]. Default 40891.5.")};
      fhicl::OptionalAtom<double> crystalRadius{ Name("crystalRadius"), Comment("LaBr crystal radius [mm]. Default 38.1.")};
      fhicl::OptionalAtom<double> crystalHalfLength{ Name("crystalHalfLength"), Comment("LaBr crystal half-length [mm]. Default 38.1.")};

      // Placement: either a distance along the crystal axis, or an explicit position.
      fhicl::OptionalAtom<double> distanceFromCrystalFace{ Name("distanceFromCrystalFace"), Comment("Place the source on-axis, this far upstream of the crystal's front face [mm]. Overrides x/y/z.")};
      fhicl::OptionalAtom<std::string> coordinates{ Name("coordinates"), Comment("Frame of x/y/z: 'global' (Mu2e, default) or 'crystal' (offset from the crystal centre).")};
      fhicl::OptionalAtom<double> x{ Name("x"), Comment("Explicit source x [mm]")};
      fhicl::OptionalAtom<double> y{ Name("y"), Comment("Explicit source y [mm]")};
      fhicl::OptionalAtom<double> z{ Name("z"), Comment("Explicit source z [mm]")};

      fhicl::OptionalAtom<std::string> emissionMode{ Name("emissionMode"), Comment("'cone' (default, solid-angle biased toward the crystal) or 'isotropic' (unbiased 4pi)")};
      fhicl::OptionalAtom<double> coneSafetyFactor{ Name("coneSafetyFactor"), Comment("Enlarges the bias cone beyond the crystal's bounding sphere. Default 1.25.")};
      fhicl::OptionalAtom<int> verbosityLevel{ Name("verbosityLevel"), Comment("Verbosity. Default 0.")};
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit CalibrationSource(const Parameters& conf);
    void produce(art::Event& event) override;
  private:
    enum Isotope { kCs137, kCo60, kNa22, kMn54, kCustom };

    CLHEP::Hep3Vector isotropicDir();
    CLHEP::Hep3Vector coneDir();
    bool acceptDirection();                       // true if this gamma is kept under biasing
    void emitGamma(double E_MeV, const CLHEP::Hep3Vector& pos, double t, GenParticleCollection& out);
    void emitDecay(const CLHEP::Hep3Vector& pos, double t, GenParticleCollection& out);

    Isotope _isotope = kCs137;
    std::string _isotopeName;
    double _windowTime = 0.0;
    double _meanDecaysPerEvent = 0.0;
    CLHEP::Hep3Vector _pos;                       // source position
    CLHEP::Hep3Vector _coneAxis;                  // unit vector, source -> crystal centre
    bool _useCone = true;
    double _cosThetaMax = -1.0;                   // cone edge; -1 means full 4pi
    double _acceptFraction = 1.0;                 // f = (1 - cosThetaMax)/2
    std::vector<double> _customE, _customProb;
    int _verbosityLevel = 0;

    art::RandomNumberGenerator::base_engine_t& _engine;
    CLHEP::RandFlat _flat;
    CLHEP::RandPoissonQ _poisson;
  };

  CalibrationSource::CalibrationSource(const Parameters& conf):
    art::EDProducer(conf),
    _windowTime(conf().windowTime()),
    _engine(createEngine(art::ServiceHandle<SeedService>()->getSeed())),
    _flat(_engine),
    _poisson(_engine, 1.0) {
      produces<GenParticleCollection>();

      _isotopeName = conf().isotope();
      if      (_isotopeName == "Cs137")  _isotope = kCs137;
      else if (_isotopeName == "Co60")   _isotope = kCo60;
      else if (_isotopeName == "Na22")   _isotope = kNa22;
      else if (_isotopeName == "Mn54")   _isotope = kMn54;
      else if (_isotopeName == "custom") _isotope = kCustom;
      else
        throw cet::exception("RANGE") << "unknown isotope '" << _isotopeName
          << "'; use Cs137, Co60, Na22, Mn54 or custom\n";

      if (_isotope == kCustom) {
        if (!conf().lineEnergies(_customE) || _customE.empty())
          throw cet::exception("RANGE") << "isotope 'custom' requires lineEnergies\n";
        if (!conf().lineIntensities(_customProb) || _customProb.size() != _customE.size())
          throw cet::exception("RANGE") << "lineIntensities must be given and match lineEnergies in length\n";
        for (size_t i = 0; i < _customE.size(); ++i) {
          if (_customE[i] <= 0.0)
            throw cet::exception("RANGE") << "custom line energies must be > 0\n";
          if (_customProb[i] < 0.0 || _customProb[i] > 1.0)
            throw cet::exception("RANGE") << "lineIntensities are per-decay probabilities and must be in [0,1]\n";
        }
      }

      double activityScale = conf().activityScale() ? *(conf().activityScale()) : 1.0;
      _verbosityLevel      = conf().verbosityLevel() ? *(conf().verbosityLevel()) : 0;
      if (_windowTime <= 0.0)
        throw cet::exception("RANGE") << "windowTime must be > 0\n";
      if (activityScale < 0.0)
        throw cet::exception("RANGE") << "activityScale must be >= 0\n";

      // Activity at the reference date, in Bq.
      double refActivity_Bq = 0.0;
      if (conf().activity() && conf().activity_uCi())
        throw cet::exception("RANGE") << "give activity OR activity_uCi, not both\n";
      if (conf().activity())          refActivity_Bq = *(conf().activity());
      else if (conf().activity_uCi()) refActivity_Bq = *(conf().activity_uCi()) * 37000.0;
      else throw cet::exception("RANGE") << "set activity [Bq] or activity_uCi\n";
      if (refActivity_Bq < 0.0)
        throw cet::exception("RANGE") << "activity must be >= 0\n";

      // Decay correction from the reference date to the measurement.
      const double elapsedDays = conf().elapsedDays() ? *(conf().elapsedDays()) : 0.0;
      double halfLife_days = 0.0;
      switch (_isotope) {
        case kCs137: halfLife_days = 30.08 * 365.25;  break;
        case kCo60:  halfLife_days = 1925.28;         break;
        case kNa22:  halfLife_days = 2.6018 * 365.25; break;
        case kMn54:  halfLife_days = 312.20;          break;
        case kCustom: break;
      }
      if (conf().halfLife_days()) halfLife_days = *(conf().halfLife_days());
      double decayFactor = 1.0;
      if (elapsedDays != 0.0) {
        if (halfLife_days <= 0.0)
          throw cet::exception("RANGE") << "elapsedDays needs halfLife_days for this isotope\n";
        decayFactor = std::exp(-std::log(2.0) * elapsedDays / halfLife_days);
      }

      double activity_Bq = refActivity_Bq * decayFactor * activityScale;
      _meanDecaysPerEvent = activity_Bq * _windowTime * 1e-9;

      // Crystal, for placement and for sizing the bias cone.
      double cX = conf().crystalX()          ? *(conf().crystalX())          : -3863.36;
      double cY = conf().crystalY()          ? *(conf().crystalY())          : 0.0;
      double cZ = conf().crystalZ()          ? *(conf().crystalZ())          : 40891.50;
      double cR = conf().crystalRadius()     ? *(conf().crystalRadius())     : 38.1;
      double cH = conf().crystalHalfLength() ? *(conf().crystalHalfLength()) : 38.1;
      if (cR <= 0.0 || cH <= 0.0)
        throw cet::exception("RANGE") << "crystalRadius and crystalHalfLength must be > 0\n";
      const CLHEP::Hep3Vector crystalCentre(cX, cY, cZ);

      // Placement. distanceFromCrystalFace puts the source on the crystal axis, upstream
      // of the front face, which is the usual bench-calibration arrangement.
      double dFace = 0.0;
      if (conf().distanceFromCrystalFace(dFace)) {
        if (dFace <= 0.0)
          throw cet::exception("RANGE") << "distanceFromCrystalFace must be > 0\n";
        _pos.set(cX, cY, cZ - cH - dFace);
      } else {
        double px = 0.0, py = 0.0, pz = 0.0;
        if (!conf().x(px) || !conf().y(py) || !conf().z(pz))
          throw cet::exception("RANGE") << "set either distanceFromCrystalFace or all of x, y, z\n";
        const std::string frame = conf().coordinates() ? *(conf().coordinates()) : "global";
        if (frame == "global")
          _pos.set(px, py, pz);
        else if (frame == "crystal")
          _pos = crystalCentre + CLHEP::Hep3Vector(px, py, pz);   // crystal has identity rotation
        else
          throw cet::exception("RANGE") << "coordinates must be 'global' or 'crystal'\n";
        const double dx = _pos.x() - cX, dy = _pos.y() - cY;
        if (std::abs(_pos.z() - cZ) <= cH && std::sqrt(dx * dx + dy * dy) <= cR)
          throw cet::exception("RANGE") << "source position " << _pos << " is inside the crystal\n";
      }

      // Bias cone: a sphere bounding the whole crystal, seen from the source. Using the
      // bounding sphere rather than the front face keeps this correct for off-axis
      // placements and for gammas entering through the side.
      const CLHEP::Hep3Vector toCrystal = crystalCentre - _pos;
      const double d = toCrystal.mag();
      const double rBound = std::sqrt(cR * cR + cH * cH);
      const double safety = conf().coneSafetyFactor() ? *(conf().coneSafetyFactor()) : 1.25;
      if (safety < 1.0)
        throw cet::exception("RANGE") << "coneSafetyFactor must be >= 1\n";

      std::string mode = conf().emissionMode() ? *(conf().emissionMode()) : "cone";
      if (mode != "cone" && mode != "isotropic")
        throw cet::exception("RANGE") << "emissionMode must be 'cone' or 'isotropic'\n";
      _useCone = (mode == "cone");

      if (d <= 0.0)
        throw cet::exception("RANGE") << "source position coincides with the crystal centre\n";
      _coneAxis = toCrystal.unit();

      const double sinThetaMax = safety * rBound / d;
      if (!_useCone || sinThetaMax >= 1.0) {
        // Source is inside or very close to the crystal: the cone would be all of 4pi,
        // so fall back to unbiased emission rather than pretend to save anything.
        _useCone = false;
        _cosThetaMax = -1.0;
        _acceptFraction = 1.0;
      } else {
        _cosThetaMax = std::sqrt(1.0 - sinThetaMax * sinThetaMax);
        _acceptFraction = 0.5 * (1.0 - _cosThetaMax);
      }

      if (_verbosityLevel > 0) {
        // On-axis disc solid angle, reported as an independent analytic cross-check of
        // the rate the simulation should produce. Only meaningful for an on-axis source
        // upstream of the front face.
        const double dFaceActual = (cZ - cH) - _pos.z();
        const bool onAxis = (std::abs(_pos.x() - cX) < 1e-6) && (std::abs(_pos.y() - cY) < 1e-6)
                            && (dFaceActual > 0.0);
        std::ostringstream extra;
        if (onAxis) {
          const double omegaDisc = 0.5 * (1.0 - dFaceActual / std::sqrt(dFaceActual * dFaceActual + cR * cR));
          extra << "\n  front-face solid angle Omega/4pi = " << omegaDisc
                << "  (analytic, on-axis disc of radius " << cR << " mm at " << dFaceActual << " mm)"
                << "\n  gammas/s striking the front face  = " << activity_Bq * omegaDisc
                << "  (per unit gamma yield; multiply by the line's per-decay probability)";
        }
        mf::LogInfo("CalibrationSource")
          << "isotope = " << _isotopeName << ", activity = " << activity_Bq << " Bq"
          << (activityScale != 1.0 ? " (scaled)" : "") << "\n"
          << "  reference activity = " << refActivity_Bq << " Bq, elapsed " << elapsedDays
          << " d, decay factor = " << decayFactor << "\n"
          << "  source position = " << _pos << " mm (crystal-relative "
          << (_pos - crystalCentre) << ")\n"
          << "  distance to crystal centre = " << d << " mm\n"
          << "  mean decays per " << _windowTime << " ns event = " << _meanDecaysPerEvent << "\n"
          << "  emission mode = " << (_useCone ? "cone (biased)" : "isotropic (unbiased)")
          << (_useCone ? "" : " -- every gamma simulated") << "\n"
          << "  bias cone half-angle = " << (_useCone ? std::acos(_cosThetaMax) * 180.0 / M_PI : 180.0) << " deg"
          << ", acceptance fraction f = " << _acceptFraction
          << (_useCone ? "  (CPU saving ~1/f)" : "")
          << extra.str();
      }
    };

  CLHEP::Hep3Vector CalibrationSource::isotropicDir() {
    double cosTheta = 2.0 * _flat.fire() - 1.0;
    double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
    double phi = 2.0 * M_PI * _flat.fire();
    return CLHEP::Hep3Vector(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
  }

  CLHEP::Hep3Vector CalibrationSource::coneDir() {
    // Uniform inside a cone of half-angle thetaMax about _coneAxis: cos(theta) uniform in
    // [cosThetaMax, 1] reproduces the correct solid-angle density.
    double cosTheta = _cosThetaMax + (1.0 - _cosThetaMax) * _flat.fire();
    double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
    double phi = 2.0 * M_PI * _flat.fire();
    CLHEP::Hep3Vector dir(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
    dir.rotateUz(_coneAxis);   // rotate the +z-referenced cone onto the real axis
    return dir;
  }

  bool CalibrationSource::acceptDirection() {
    // Under biasing, a gamma survives with probability f and is then thrown inside the
    // cone. Applying this INDEPENDENTLY per gamma keeps cascade coincidences right: for
    // Co-60 both gammas must be accepted for a sum-peak event, probability f^2, which is
    // what isotropic emission of two independent directions would also give.
    return !_useCone || (_flat.fire() < _acceptFraction);
  }

  void CalibrationSource::emitGamma(double E_MeV, const CLHEP::Hep3Vector& pos, double t,
                                    GenParticleCollection& out) {
    if (!acceptDirection())
      return;
    CLHEP::Hep3Vector dir = _useCone ? coneDir() : isotropicDir();
    CLHEP::HepLorentzVector mom(E_MeV * dir, E_MeV);   // photon: |p| = E
    out.push_back(GenParticle(PDGCode::gamma, GenId::particleGun, pos, mom, t));
  }

  void CalibrationSource::emitDecay(const CLHEP::Hep3Vector& pos, double t,
                                    GenParticleCollection& out) {
    switch (_isotope) {
      case kCs137:
        if (_flat.fire() < kCs137_GammaProb) emitGamma(kCs137_Gamma, pos, t, out);
        if (_flat.fire() < kCs137_KaProb)    emitGamma(kCs137_KaE,   pos, t, out);
        if (_flat.fire() < kCs137_KbProb)    emitGamma(kCs137_KbE,   pos, t, out);
        break;

      case kCo60:
        // Same vertex and time: a true cascade, so both can be absorbed together and
        // produce the 2505 keV sum peak.
        if (_flat.fire() < kCo60_Gamma1Prob) emitGamma(kCo60_Gamma1, pos, t, out);
        if (_flat.fire() < kCo60_Gamma2Prob) emitGamma(kCo60_Gamma2, pos, t, out);
        break;

      case kNa22:
        if (_flat.fire() < kNa22_GammaProb) emitGamma(kNa22_Gamma, pos, t, out);
        if (_flat.fire() < kNa22_BetaPlusBr) {
          // Two annihilation quanta. They are physically back to back, but a single
          // detector can only ever intercept one of them, so emitting two independent
          // isotropic gammas gives the identical singles rate. (This approximation would
          // NOT be valid for a two-detector coincidence study.)
          emitGamma(kNa22_Annih, pos, t, out);
          emitGamma(kNa22_Annih, pos, t, out);
        }
        break;

      case kMn54:
        if (_flat.fire() < kMn54_GammaProb) emitGamma(kMn54_Gamma, pos, t, out);
        break;

      case kCustom:
        for (size_t i = 0; i < _customE.size(); ++i)
          if (_flat.fire() < _customProb[i]) emitGamma(_customE[i], pos, t, out);
        break;
    }
  }

  void CalibrationSource::produce(art::Event& event) {
    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);
    int nDecays = _poisson.fire(_meanDecaysPerEvent);
    for (int i = 0; i < nDecays; ++i) {
      double t = _flat.fire() * _windowTime;   // uniform in window: constant activity
      emitDecay(_pos, t, *output);
    }
    event.put(std::move(output));
  };
}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::CalibrationSource)
