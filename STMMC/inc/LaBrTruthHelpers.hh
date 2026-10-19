#ifndef STMMC_LaBrTruthHelpers_hh
#define STMMC_LaBrTruthHelpers_hh
//
// Helpers shared by the LaBr digitizer and the LaBr truth modules, so that all three agree
// on (1) how much scintillation light a step makes and (2) which source a primary came from.
//
// 1) Alpha light quenching. LaBr3 produces ~0.3x less light per MeV for alphas than for
//    electrons, so the 5-7.4 MeV Ac-227 chain alphas appear at ~1.5-2.4 MeV electron-equivalent.
//    Every step made by an alpha is scaled by
//        QF(E0) = alphaQuenchingFactor + alphaQuenchingSlope * (E0 - 6 MeV)
//    where E0 is the alpha's initial kinetic energy (the light output depends on the whole
//    track, not the local step). With the defaults (1, 0) nothing is changed. Fit both to the
//    positions of the alpha peaks in bench data.
//
// 2) Source classification of a primary:
//      GenId::STMAc227                      -> STMMCPulse::Ac227
//      GenId::STMLa138                      -> STMMCPulse::La138
//      otherwise, vertex inside the crystal -> STMMCPulse::La138  (files made before the tags)
//      otherwise                            -> STMMCPulse::photonGun (any EXTERNAL source:
//                                              photon gun, calibration source, beam)
//    This replaces the old start-z threshold, which failed for a source placed beside the
//    crystal rather than upstream of it.
//
// Original author: Sena Clarke
//

#include <algorithm>
#include <cmath>
#include <tuple>

#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/STMMCPulse.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"

namespace mu2e {
  namespace LaBrTruth {

    constexpr double kQuenchingReferenceEnergy = 6.0;   // MeV

    // Light-equivalent energy of a step [same units as ionizingEdep].
    inline double lightEquivalentEdep(const StepPointMC& step, double alphaQF, double alphaSlope) {
      const double edep = step.ionizingEdep();
      if (alphaQF == 1.0 && alphaSlope == 0.0)
        return edep;   // do not touch the SimParticle Ptr unless quenching is enabled
      const SimParticle* sp = step.simParticle().get();
      if (!sp || sp->pdgId() != PDGCode::alpha)
        return edep;
      const CLHEP::HepLorentzVector& p4 = sp->startMomentum();
      const double E0 = p4.e() - p4.m();
      const double qf = alphaQF + alphaSlope * (E0 - kQuenchingReferenceEnergy);
      return edep * std::max(0.0, qf);
    }

    // Walk the parent chain to the true primary. SimParticle::originParticle() only recurses
    // for the self-parent special case, so it returns ordinary secondaries themselves.
    inline const SimParticle* findPrimary(const StepPointMC& step) {
      const SimParticle* primary = step.simParticle().get();
      while (primary && primary->hasParent()) {
        const SimParticle* up = primary->parent().get();
        if (!up) break;   // truncated chain: keep the lowest available ancestor
        primary = up;
      }
      return primary;
    }

    struct CrystalCylinder {
      double x, y, z;        // centre [mm]
      double radius;         // [mm]
      double halfLength;     // [mm]
      double tolerance = 0.1;
      bool contains(const CLHEP::Hep3Vector& p) const {
        const double dx = p.x() - x, dy = p.y() - y;
        return std::abs(p.z() - z) <= halfLength + tolerance &&
               std::sqrt(dx * dx + dy * dy) <= radius + tolerance;
      }
    };

    // Identifies one decay. All products of a decay are generated at the same vertex and time,
    // so primaries sharing (source, start vertex, start time) are one decay and form ONE pulse
    // in the detector. Rounded to 1 um and 1 ps. A photon-gun photon is its own "decay".
    using DecayKey = std::tuple<int, long long, long long, long long, long long>;
    inline DecayKey decayKey(const SimParticle& primary, int source) {
      const CLHEP::Hep3Vector v = primary.startPosition();
      return DecayKey{source, std::llround(v.x() * 1e3), std::llround(v.y() * 1e3),
                      std::llround(v.z() * 1e3), std::llround(primary.startGlobalTime() * 1e3)};
    }

    inline int classifyPrimary(const SimParticle& primary, const CrystalCylinder& crystal) {
      const auto& gp = primary.genParticle();
      if (gp.isNonnull() && gp.isAvailable()) {
        const GenId id = gp->generatorId();
        if (id == GenId::STMAc227) return STMMCPulse::Ac227;
        if (id == GenId::STMLa138) return STMMCPulse::La138;
      }
      return crystal.contains(primary.startPosition()) ? STMMCPulse::La138
                                                       : STMMCPulse::photonGun;
    }

  } // namespace LaBrTruth
} // namespace mu2e
#endif
