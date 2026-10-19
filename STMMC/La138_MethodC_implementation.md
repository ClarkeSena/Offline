# Implementing La-138 Intrinsic Activity via Direct Decay-Product Generation (Method C)

> **STATUS (2026-06-17): IMPLEMENTED.** This was the design/implementation guide; Method C
> is now in `La138Source_module.cc` and standalone-validated. For current overall state,
> the file map, what's verified vs open, and gotchas, read **`La138_STATUS.md`** first.
> This document remains the reference for *why* Method C and the decay-scheme physics.

**Audience:** a fresh Claude Code session (or developer) with no prior context on this task.
**Goal:** make the LaBr3 simulation reproduce the intrinsic ¹³⁸La self-scintillation
peaks that appear in real data but not in simulation, by having the generator emit
the decay products directly (Method C), bypassing Geant4 radioactive decay.

This document is self-contained: it explains *why* this approach, the *physics*, the
*exact code* to write, how to *build/run*, how to *validate*, and every *pitfall*
discovered while scoping it.

---

## 0. Repository orientation (read first)

- Working dir: `/Users/senaclarke/Desktop/Physics/Mu2e/Simulation` (edit here, then copy
  to the gpvm at `/exp/mu2e/app/users/sclarke/Simulation/` to build/run — **the build does
  not happen locally**).
- All STM code: `Offline/STMMC/`. Modules in `src/`, fcl in `fcl/`, build in `CMakeLists.txt`.
- Build: `muse setup` then `muse build -j8`. Run: `mu2e -c <fcl>`.
- **Build pitfalls (verified, must obey):**
  - Module source files must be named `*_module.cc`.
  - `SeedService` linking fails for STMMC modules — use `std::default_random_engine`,
    NOT `art::ServiceHandle<SeedService>`. (The existing module already does this.)
  - FHiCL files must be pure ASCII — no em dashes, smart quotes, non-breaking spaces.
- Mu2e unit conventions in `GenParticle`: **energy/momentum in MeV, position in mm,
  time in ns**, all in **Mu2e coordinates** (the primary generator adds `mu2eOrigin`).

---

## 1. Why Method C (background context)

The simulation only ever fires *external* particles (photon gun, beam backgrounds), so
nothing radioactive happens *inside* the crystal — hence no intrinsic peaks. ¹³⁸La is
0.09% of the lanthanum in the LaBr₃ lattice and decays, producing characteristic peaks.

Three ways to simulate it were considered:

- **A — inject La-138 ion + analogue `G4RadioactiveDecay`.** *Fails.* La-138's half-life
  is 1.02×10¹¹ yr, so analogue Monte Carlo samples a decay time ~10²⁷ ns; the in-window
  (1695 ns) decay probability is ~3.7×10⁻²⁵. The ion never decays in-event → zero signal.
- **B — inject ion + biased decay (`G4Radioactivation` / `Mu2eG4BiasedRDPhysics`,
  `radiationVRmode`).** Works for an energy *spectrum*, but every decay carries a
  statistical weight ~10⁻²⁵. Weights work in a histogram, **not** in a continuous
  waveform/time-series — you can't put "10⁻²⁵ of a pulse" in the ADC stream at the right
  rate. So B is mismatched to this waveform pipeline.
- **C — generator emits the decay products directly (THIS DOCUMENT).** Unit weight,
  correct rate, drops straight into the waveform. You hand-code only the *source*
  emission; Geant4 still does all *transport* (Compton, photoabsorption, escape, and the
  coincidence summing in the crystal), so the photopeak/Compton/escape/sum-peak structure
  emerges correctly. You own only: branching ratios, γ line energies, β spectrum shape,
  X-ray energy/yield.

---

## 2. The ¹³⁸La decay scheme (the physics you are coding)

VERIFY these against ENSDF/NNDC before publication; values below are good working numbers.

```
                 138La  (T1/2 = 1.02e11 yr)
                  /                         \
   EC  ~66.4%    /                           \   beta-  ~33.6%
               v                              v
            138Ba* (1435.8 keV)            138Ce* (788.7 keV)
               |                              |
   gamma 1435.795 keV                gamma 788.742 keV
   + Ba K X-ray (~32-37 keV)         + beta- continuum (endpoint ~263 keV)
   + neutrino (INVISIBLE)            + antineutrino (INVISIBLE)
```

What deposits in the crystal (everything except the neutrino is absorbed locally because
the crystal is large compared to the ranges):

- **EC branch:** 1436 keV γ + a Ba K X-ray. Both absorbed -> **sum peak ~1468 keV**;
  if the ~32 keV X-ray escapes (near a surface) -> **1436 keV photopeak**. Neutrino escapes.
- **beta- branch:** 789 keV γ + the β electron (range << crystal, fully absorbed) ->
  a **789 keV peak riding on a continuum** up to ~789+263 ~ 1052 keV.

Key numbers (MeV unless noted):

| Quantity | Value | Notes |
|---|---|---|
| EC branch fraction | 0.664 | rest (0.336) is β⁻ |
| EC γ energy | 1.435795 | ¹³⁸Ba* → ¹³⁸Ba |
| β⁻ γ energy | 0.788742 | ¹³⁸Ce* → ¹³⁸Ce |
| β⁻ endpoint (kinetic) | 0.263 | = Q_β⁻(~1052 keV) − 789 keV |
| electron mass | 0.510999 | |
| Ba Kα X-ray | 0.0322 | ~83% of K X-rays |
| Ba Kβ X-ray | 0.0364 | ~17% of K X-rays |
| P(K X-ray per EC) | 0.71 | = P_K(~0.80) × ω_K(Ba, ~0.894); else Auger (still local) |

β⁻ spectrum (allowed shape): `N(E) dE ∝ p E (E0 − E)^2 dE`, with E = total energy =
KE + m_e, p = √(E²−m_e²), E0 = endpoint total energy. The Fermi function F(Z=58, E)
(Z of the Ce daughter) is an accuracy refinement, omitted in the first pass.

---

## 3. Current state of the files

These already exist from the prior session:

- `Offline/STMMC/src/La138Source_module.cc` — **currently injects La-138 IONS** (relies
  on radioactive decay). This is what you REWRITE for Method C. It already has the correct
  activity calculation (computes ~544 Bq from crystal geometry/material) and uniform-in-
  cylinder vertex sampling — KEEP those. You replace only the ion emission with product
  emission, and remove the now-unused `La138PDG` constant.
- `Offline/STMMC/src/GenParticleMerger_module.cc` — concatenates GenParticleCollections so
  PhotonGun + La138Source feed one g4run. KEEP as-is.
- `Offline/STMMC/CMakeLists.txt` — already registers `La138Source` and `GenParticleMerger`
  with `Offline::MCDataProducts` and `Offline::DataProducts`. NO CMake change needed.
- `Offline/STMMC/fcl/PhotonGunPlusLa138ToLaBrWaveform.fcl` — combined pipeline. Currently
  sets `turnOnRadioactiveDecay: true`; for Method C this must be **false** (see step 5).

**Geometry facts (use these; do NOT query GeometryService):** the LaBr crystal is a solid
cylinder, CrystalR = 38.1 mm, CrystalL = 76.2 mm, centre in Mu2e coords
(−3863.36, 0, 40891.50) mm. NOTE: `STMGeometryDump` / `labr->originInMu2e()` report the
crystal z as 40993.1 — that is **101.6 mm (4 in, the front tungsten shield) too high** due
to a reference-anchor mismatch between `STMMaker.cc` and `constructSTM.cc`. The value
40891.50 matches the actual G4 step positions and is what the digitizer uses. Querying the
geometry service would place decays 101.6 mm off (partly outside the crystal), so the
hard-coded 40891.50 is correct.

---

## 4. Rewrite `La138Source_module.cc` (the core step)

Replace the ENTIRE file with the following. It keeps the constructor's activity
calculation and vertex sampling, adds decay-scheme constants + helper methods, and
rewrites `produce()` to emit branched products. (Diff-conscious: the constructor body is
unchanged except for removing the `La138PDG` usage; the class declaration gains helper
method prototypes and a `betaShapeMax_` member.)

```cpp
// Intrinsic LaBr3 radioactivity source (Method C): for each La-138 decay (drawn at the
// rate implied by the crystal's lanthanum content), emit the decay PRODUCTS directly as
// GenParticles at a uniform vertex inside the crystal. Geant4 then transports them, so
// the photopeak/Compton/escape structure and the EC coincidence sum (1436 + Ba X-ray ~
// 1468 keV) and the beta- 789 keV-plus-continuum emerge from transport. No G4 radioactive
// decay is used (turnOnRadioactiveDecay must be false), so events are unit-weight and land
// in the waveform at the correct rate. Merge with PhotonGun via GenParticleMerger.
// Original author: Sena Clarke

#include <cmath>
#include <random>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"

#include "cetlib_except/exception.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"

namespace {
  // ---- 138La decay scheme (ENSDF/NNDC). VERIFY/REFINE before publication. ----
  // Energies in MeV.
  constexpr double kEcFraction       = 0.664;      // EC branch fraction (rest is beta-)
  constexpr double kGammaEC_MeV      = 1.435795;   // 138Ba* -> 138Ba gamma (EC branch)
  constexpr double kGammaBeta_MeV    = 0.788742;   // 138Ce* -> 138Ce gamma (beta- branch)
  constexpr double kBetaEndpoint_MeV = 0.263;      // beta- kinetic endpoint = Qb- - 789 keV
  constexpr double kElectronMass_MeV = 0.510999;
  constexpr double kXrayProb         = 0.71;       // P(K X-ray | EC) = P_K * omega_K(Ba)
  constexpr double kKalpha_MeV       = 0.0322;     // Ba K-alpha (weighted)
  constexpr double kKbeta_MeV        = 0.0364;     // Ba K-beta
  constexpr double kKalphaFrac       = 0.83;       // fraction of K X-rays that are K-alpha
}

namespace mu2e {
  class La138Source : public art::EDProducer {
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config {
      fhicl::Atom<double> x{ Name("x"), Comment("Crystal centre x [mm]")};
      fhicl::Atom<double> y{ Name("y"), Comment("Crystal centre y [mm]")};
      fhicl::Atom<double> z{ Name("z"), Comment("Crystal centre z [mm]")};
      fhicl::Atom<double> radius{ Name("radius"), Comment("Crystal radius [mm]")};
      fhicl::Atom<double> halfLength{ Name("halfLength"), Comment("Crystal half-length along z [mm]")};
      fhicl::OptionalAtom<double> density{ Name("density"), Comment("Crystal density [g/cm^3]. Default 5.08 (LaBr3).")};
      fhicl::OptionalAtom<double> molarMass{ Name("molarMass"), Comment("Molar mass [g/mol]. Default 378.6 (LaBr3).")};
      fhicl::OptionalAtom<double> laPerMolecule{ Name("laPerMolecule"), Comment("La atoms per formula unit. Default 1.")};
      fhicl::OptionalAtom<double> la138Abundance{ Name("la138Abundance"), Comment("La-138 isotopic abundance. Default 0.0009.")};
      fhicl::OptionalAtom<double> halfLife_yr{ Name("halfLife_yr"), Comment("La-138 half-life [yr]. Default 1.02e11.")};
      fhicl::Atom<double> windowTime{ Name("windowTime"), Comment("Event time window [ns]; match micropulseTime.")};
      fhicl::OptionalAtom<double> activityScale{ Name("activityScale"), Comment("Multiplier on physical activity for statistics. Default 1.")};
      fhicl::OptionalAtom<int> verbosityLevel{ Name("verbosityLevel"), Comment("Verbosity. Default 0.")};
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit La138Source(const Parameters& conf);
    void produce(art::Event& event) override;
  private:
    CLHEP::Hep3Vector isotropicDir();
    void emitGamma(double E_MeV, const CLHEP::Hep3Vector& pos, double t, GenParticleCollection& out);
    void emitElectron(double KE_MeV, const CLHEP::Hep3Vector& pos, double t, GenParticleCollection& out);
    double sampleBetaKE();
    double sampleBaKXray();

    double x = 0, y = 0, z = 0, radius = 0, halfLength = 0;
    double windowTime = 0;
    double meanDecaysPerEvent = 0;
    double betaShapeMax_ = 0;   // precomputed rejection-sampling bound for the beta shape
    int verbosityLevel = 0;

    std::default_random_engine _randomGen;   // SeedService linking workaround (STMMC)
    std::poisson_distribution<int> _poisson;
    std::uniform_real_distribution<double> _uniform{0.0, 1.0};
  };

  La138Source::La138Source(const Parameters& conf):
    art::EDProducer(conf),
    x(conf().x()), y(conf().y()), z(conf().z()),
    radius(conf().radius()), halfLength(conf().halfLength()),
    windowTime(conf().windowTime()) {
      produces<GenParticleCollection>();

      double density       = conf().density()        ? *(conf().density())        : 5.08;
      double molarMass     = conf().molarMass()      ? *(conf().molarMass())      : 378.6;
      double laPerMolecule = conf().laPerMolecule()  ? *(conf().laPerMolecule())  : 1.0;
      double la138Abund    = conf().la138Abundance() ? *(conf().la138Abundance()) : 0.0009;
      double halfLife_yr   = conf().halfLife_yr()    ? *(conf().halfLife_yr())    : 1.02e11;
      double activityScale = conf().activityScale()  ? *(conf().activityScale())  : 1.0;
      verbosityLevel       = conf().verbosityLevel() ? *(conf().verbosityLevel()) : 0;

      if (radius <= 0 || halfLength <= 0 || windowTime <= 0)
        throw cet::exception("RANGE") << "radius, halfLength and windowTime must be > 0\n";

      // Activity from the lanthanum content of the crystal volume.
      const double NA = 6.02214076e23;
      const double secPerYear = 3.1557600e7;
      double V_cm3 = M_PI * std::pow(radius / 10.0, 2) * (2.0 * halfLength / 10.0); // mm->cm
      double mass_g = density * V_cm3;
      double nLa138 = (mass_g / molarMass) * NA * laPerMolecule * la138Abund;
      double lambda = std::log(2.0) / (halfLife_yr * secPerYear);
      double activity_Bq = lambda * nLa138 * activityScale;
      meanDecaysPerEvent = activity_Bq * windowTime * 1e-9;
      _poisson = std::poisson_distribution<int>(meanDecaysPerEvent);

      // Precompute the beta-shape maximum for rejection sampling.
      for (int i = 1; i < 1000; ++i) {
        double KE = kBetaEndpoint_MeV * i / 1000.0;
        double E = KE + kElectronMass_MeV;
        double p = std::sqrt(E * E - kElectronMass_MeV * kElectronMass_MeV);
        double E0 = kBetaEndpoint_MeV + kElectronMass_MeV;
        double s = p * E * std::pow(E0 - E, 2);
        if (s > betaShapeMax_) betaShapeMax_ = s;
      }

      if (verbosityLevel > 0) {
        mf::LogInfo("La138Source")
          << "crystal V = " << V_cm3 << " cm^3, mass = " << mass_g << " g\n"
          << "N(La-138) = " << nLa138 << ", activity = " << activity_Bq << " Bq"
          << (activityScale != 1.0 ? " (scaled)" : "") << "\n"
          << "mean decays per " << windowTime << " ns event = " << meanDecaysPerEvent;
      }
    };

  CLHEP::Hep3Vector La138Source::isotropicDir() {
    double cosTheta = 2.0 * _uniform(_randomGen) - 1.0;
    double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
    double phi = 2.0 * M_PI * _uniform(_randomGen);
    return CLHEP::Hep3Vector(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta);
  }

  void La138Source::emitGamma(double E_MeV, const CLHEP::Hep3Vector& pos, double t,
                              GenParticleCollection& out) {
    CLHEP::Hep3Vector dir = isotropicDir();
    CLHEP::HepLorentzVector mom(E_MeV * dir, E_MeV);   // photon: |p| = E
    out.push_back(GenParticle(PDGCode::gamma, GenId::particleGun, pos, mom, t));
  }

  void La138Source::emitElectron(double KE_MeV, const CLHEP::Hep3Vector& pos, double t,
                                 GenParticleCollection& out) {
    double E = KE_MeV + kElectronMass_MeV;
    double p = std::sqrt(std::max(0.0, E * E - kElectronMass_MeV * kElectronMass_MeV));
    CLHEP::Hep3Vector dir = isotropicDir();
    CLHEP::HepLorentzVector mom(p * dir, E);
    out.push_back(GenParticle(PDGCode::e_minus, GenId::particleGun, pos, mom, t));
  }

  double La138Source::sampleBetaKE() {
    const double m = kElectronMass_MeV;
    const double E0 = kBetaEndpoint_MeV + m;
    while (true) {
      double KE = _uniform(_randomGen) * kBetaEndpoint_MeV;
      double E = KE + m;
      double p = std::sqrt(E * E - m * m);
      double s = p * E * std::pow(E0 - E, 2);
      if (_uniform(_randomGen) * betaShapeMax_ <= s) return KE;
    }
  }

  double La138Source::sampleBaKXray() {
    return (_uniform(_randomGen) < kKalphaFrac) ? kKalpha_MeV : kKbeta_MeV;
  }

  void La138Source::produce(art::Event& event) {
    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);
    int nDecays = _poisson(_randomGen);
    for (int i = 0; i < nDecays; ++i) {
      // Uniform vertex in the homogeneous crystal (axis along z). Because the crystal is
      // a homogeneous LaBr3 cylinder, uniform-per-volume IS La-density-weighted; the
      // sqrt() makes the radial sampling uniform per unit area (not per unit radius).
      double r = radius * std::sqrt(_uniform(_randomGen));
      double phi = 2.0 * M_PI * _uniform(_randomGen);
      double dz = (2.0 * _uniform(_randomGen) - 1.0) * halfLength;
      CLHEP::Hep3Vector pos(x + r * std::cos(phi), y + r * std::sin(phi), z + dz);
      double t = _uniform(_randomGen) * windowTime;   // uniform in window: constant activity

      if (_uniform(_randomGen) < kEcFraction) {
        // EC branch: 1436 keV gamma (+ Ba K X-ray); neutrino is invisible.
        emitGamma(kGammaEC_MeV, pos, t, *output);
        if (_uniform(_randomGen) < kXrayProb)
          emitGamma(sampleBaKXray(), pos, t, *output);
      } else {
        // beta- branch: 789 keV gamma + beta electron (continuum).
        emitGamma(kGammaBeta_MeV, pos, t, *output);
        emitElectron(sampleBetaKE(), pos, t, *output);
      }
    }
    event.put(std::move(output));
  };
}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::La138Source)
```

**Notes on the code:**
- `HepLorentzVector(E_MeV * dir, E_MeV)` uses the `(Hep3Vector, double)` constructor;
  for a photon |p| = E so this is on-shell. For the electron, |p| = √(E²−m²).
- All products of one decay share `pos` and `t` -> they sum in the crystal (coincidence).
- `GenId::particleGun` is reused (same as PhotonGun); it is just a provenance tag.
- The neutrino is simply not emitted (it carries away energy invisibly — correct).

---

## 5. FHiCL changes

In `Offline/STMMC/fcl/PhotonGunPlusLa138ToLaBrWaveform.fcl`:

1. **Turn OFF radioactive decay** (Method C does not use it):
   ```
   physics.producers.g4run.physics.turnOnRadioactiveDecay: false
   ```
   (If the line currently says `true`, change it; or delete it — default is false.)
2. Everything else stays: the `la138` producer block (x/y/z/radius/halfLength/windowTime),
   the `GenParticleMerger` (`genCombined` with inputTags `["generate","la138"]`), and
   `physics.producers.g4run.inputs.primaryTag: "genCombined"`.
3. The crystal geometry in the `la138` block must be:
   ```
   x : -3863.36
   y : 0.0
   z : 40891.50
   radius : 38.1
   halfLength : 38.1
   windowTime : 1695.0
   ```
4. For testing/visualising the peaks quickly, raise statistics with
   `activityScale : 1000` (physical is ~544 Bq -> ~9.2e-4 decays/event, so at scale 1 you
   need ~1e6 events to get ~1000 decays). Set back to 1 (or remove) for a physically
   normalised combined spectrum.

To run La-138 ALONE (no photon gun), make a copy of the fcl, drop `generate` from the
merger inputTags and the path, and keep only `la138, genCombined, g4run, DigiLaBr`.

---

## 6. Build and run

```bash
# on the gpvm, after copying the edited files over
muse setup
muse build -j8

# combined signal + intrinsic, boosted stats for a first look
mu2e -c Offline/STMMC/fcl/PhotonGunPlusLa138ToLaBrWaveform.fcl -n 2000
```

(`CMakeLists.txt` already registers `La138Source`; no build-config change is needed.)

---

## 7. Validation (do these in order)

1. **Generator sanity:** with `verbosityLevel : 1`, confirm the `La138Source` log prints
   `activity = ~544 Bq` and `mean decays per 1695 ns event = ~9.2e-4` (at activityScale 1).
   At scale 1000 expect ~0.92 decays/event.
2. **Products are created in the crystal:** run with high `activityScale`, dump the
   GenParticleCollection or the `STMDet` StepPointMCs; confirm gammas (1.436 / 0.789 MeV),
   electrons, and ~0.032 MeV X-rays appear, and that StepPointMC vertices lie within the
   crystal cylinder (|z − 40891.5| < 38.1, r < 38.1).
3. **Truth energy spectrum** (before digitization): histogram the total `ionizingEdep`
   per decay. Expect:
   - a **1468 keV** sum peak (1436 γ + 32 keV X-ray both absorbed),
   - a **1436 keV** peak (X-ray escaped),
   - a **789 keV** peak on a **continuum up to ~1052 keV** (β summed),
   - Compton continua below each γ line.
4. **Digitized spectrum:** run the matched-filter/peak-height analysis on the waveforms;
   the same features should appear at the calibrated energies (recall the peak-height
   calibration: energy = 0.57 × ADC_peak, set via `peakKeVPerADC` in the digitizer).
5. **Rate check:** at `activityScale : 1`, the intrinsic pulse rate should be
   ~544 Hz-equivalent (9.2e-4 per 1695 ns micropulse), i.e. ~40x rarer than a 20 kHz signal.

---

## 8. Pitfalls checklist (every gotcha for this task)

- [ ] **Do NOT enable `turnOnRadioactiveDecay`** for Method C — it is unused and, for
      La-138, would not decay anything anyway (10¹¹-yr half-life; see section 1A).
- [ ] **Units are MeV** in GenParticle — 1435.795 keV = 1.435795 MeV. Getting this wrong
      by 1000x puts pulses off-scale.
- [ ] **Crystal z = 40891.50**, NOT the 40993.1 from STMGeometryDump/originInMu2e
      (101.6 mm too high). Do not "fix" this by querying GeometryService.
- [ ] **`std::default_random_engine`**, not SeedService (STMMC linking limitation).
- [ ] File stays named `La138Source_module.cc` (art plugin naming).
- [ ] FHiCL pure ASCII only.
- [ ] All products of one decay share `pos` and `t` (else no coincidence summing).
- [ ] `radius` and `halfLength` are 38.1 each (CrystalR=38.1, CrystalL=76.2 -> half=38.1).
- [ ] Remove the old `La138PDG` ion constant — Method C emits γ/e⁻, not the ion.
- [ ] The `GenParticleMerger` and `g4run.inputs.primaryTag: genCombined` must stay so the
      photon gun and La-138 products both reach Geant4.

---

## 9. Refinements / future work (optional, after the basic version validates)

- **Fermi function** F(Z=58, E) in `sampleBetaKE()` for an accurate β continuum shape.
- **Exact branching/energies** from current ENSDF (the table values are working numbers).
- **Internal conversion electrons** and **L X-rays** for fine structure.
- **²²⁷Ac contamination** (the α "bumps"): a SEPARATE, larger task. Its decay CHAIN is
  better done with biased decay (Method B) than hand-coded. Crucially, ²²⁷Ac is NOT
  uniform in the crystal (segregates during growth) — it needs a radial/axial density
  profile, not the flat sampling used for ¹³⁸La. Also requires **alpha quenching**
  (Birks' law / a particle-type-dependent light yield) in the digitizer, because Geant4
  gives the full alpha energy but a scintillator produces far less light per MeV for
  alphas — without quenching the α features land at the wrong (too-high) energy. The
  current digitizer applies one `lightYield` to all particles.
- **Absolute combined normalization:** with `activityScale : 1` the intrinsic rate is
  physical; only boost it for statistics/visualisation, then restore.

---

## 10. One-paragraph summary for the implementer

Rewrite `Offline/STMMC/src/La138Source_module.cc` (section 4) so that, per Poisson-drawn
decay, it emits the ¹³⁸La decay products as GenParticles at a uniform vertex in the
crystal cylinder: 66.4% EC -> a 1.436 MeV γ plus (71%) a ~0.032 MeV Ba X-ray; else β⁻ ->
a 0.789 MeV γ plus an electron sampled from the allowed β spectrum (endpoint 0.263 MeV);
all sharing vertex and time. Set `turnOnRadioactiveDecay: false` in
`PhotonGunPlusLa138ToLaBrWaveform.fcl`. Build, run with `activityScale` boosted, and
validate the 1468/1436/789 keV features in the truth and digitized spectra. No CMake
change. Obey the pitfalls in section 8 — especially MeV units, crystal z = 40891.50, and
`std::default_random_engine`.
