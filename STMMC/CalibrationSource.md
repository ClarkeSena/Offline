# Calibration source — validating the La-138 intrinsic rate

**Problem this solves.** `La138Source` computes the crystal's intrinsic activity from its
lanthanum content (~544 Bq) and that number is currently taken on trust, then multiplied by
an arbitrary `activityScale` for statistics. Nothing checks whether the simulation actually
*delivers* that rate into the spectrum.

**The method.** A sealed source of known activity at a known distance is a **transfer
standard**. Its rate in the detector is independently predictable from first principles —
activity × gamma yield × solid angle × intrinsic efficiency — and measurable on a bench. Run
it through the same simulation chain as La-138, at the same absolute-rate convention, and:

- if the simulation reproduces the **measured** source rate, the chain's absolute-rate
  machinery (geometry, acceptance, containment, livetime) is validated;
- the La-138 rate produced by that same validated chain can then be believed;
- and if the two disagree, the source side is known independently, so the discrepancy is
  localised rather than ambiguous.

---

## Files

| File | Role |
|---|---|
| `src/CalibrationSource_module.cc` | Generator. Emits an isotope's gammas at an absolute Bq rate. |
| `fcl/CalibrationSourceToLaBrWaveform.fcl` | Source + La-138 + G4 + digitizer + truth-pulse chain. |
| `src/analyseCalibrationRates.C` | Turns truth pulses into counts/s and prints the ratio. |

C++ changed, so this needs `muse build -j8`.

## Running it

```bash
mu2e -c Offline/STMMC/fcl/CalibrationSourceToLaBrWaveform.fcl
root -l -b -q 'Offline/STMMC/src/analyseCalibrationRates.C("nts.CalibrationSource.root",100,37000)'
```

Arguments to the macro: file, `activityScale` (**must match the fcl**), source activity in Bq,
threshold in keV (default 20), event-count fallback.

For the source alone — measuring intrinsic peak efficiency with no La-138 contamination:

```
physics.producers.la138.activityScale : 0
```

## The rate convention

Both generators use the identical convention, which is what makes the comparison meaningful:

```
mean decays per event = activity[Bq] * activityScale * windowTime[ns] * 1e-9
```

> **`calsource.activityScale` and `la138.activityScale` must be equal.** They exist only to
> buy statistics; if they differ, the ratio — the entire point of the exercise — is wrong.
> The analysis macro divides the scale back out, so everything it prints is a physical rate.

## Isotopes

| `isotope` | Lines (keV) | Notes |
|---|---|---|
| `Cs137` | 661.657 (85.1%), Ba K X-rays 32.06 (5.8%) / 36.4 (1.33%) | Default. 662 keV is where `excessNoiseFactor` is calibrated, so this run also tests the resolution model. |
| `Co60` | 1173.228 (99.85%) + 1332.492 (99.98%) | True cascade — same vertex and time, so the 2505 keV sum peak appears. Brackets the 1809 signal line. |
| `Na22` | 1274.537 (99.94%), 511 ×2 (90.3%) | Positron assumed to annihilate in the encapsulation. |
| `Mn54` | 834.848 (99.98%) | Single clean line. |
| `custom` | `lineEnergies` + `lineIntensities` | Intensities are **absolute per-decay probabilities**, not relative. |

Decay data are ENSDF/NNDC values and carry the same "verify before publication" caveat as
`La138Source`. Multi-line sources with strong cascades (Ba-133, Eu-152) can be entered as
`custom`, but their lines would be treated as **independent** — true coincidence summing is
only modelled for the built-in `Co60` cascade.

## Placement

`distanceFromCrystalFace` puts the source on-axis, that far upstream of the crystal's front
face (z = 40891.50 − 38.1 = 40853.4 mm). Or give explicit `x`, `y`, `z`.

## Solid-angle biasing

A point source emits into 4π but the crystal subtends under 1% of it at 200 mm, so 99% of the
Geant4 time would be wasted. With `emissionMode : "cone"` (the default) each gamma is
independently kept with probability *f* and, if kept, thrown uniformly inside a cone that
contains the crystal's bounding sphere.

This is mathematically identical to isotropic emission *for every gamma that could reach the
crystal*, with unit weights throughout. At 200 mm, f ≈ 0.020 — a **~49× speedup**.

Applying acceptance independently per gamma keeps cascades right: for Co-60 both gammas must
be accepted for a sum-peak event, probability *f*², exactly as two independent isotropic
directions would give.

> **What biasing discards:** a gamma emitted *away* from the detector that scatters off
> surrounding material back into it. On a bench that is a small correction; inside the full
> Mu2e geometry, with the SSC collimator and shielding nearby, it may not be. Cross-check any
> final number with `emissionMode : "isotropic"`.

## What to expect

Analytic prediction for the default configuration — 37 kBq (1 µCi) Cs-137 at 200 mm, 3″×3″
LaBr₃ (76.2 mm × 76.2 mm), µ/ρ = 0.0745 cm²/g at 662 keV:

```
  662 keV gammas emitted        31487 /s
  striking the front face         278 /s     (Omega/4pi = 0.883%)
  interacting (any deposit)       263 /s     (94.4% interact in 7.62 cm)
  in the photopeak                144 /s     (intrinsic peak eff ~0.55)

  La-138 intrinsic activity       544 Bq     (348 cm^3, 0.09% abundance)

  EXPECTED source / La-138 detected-rate ratio  ~0.48
```

That is a good design point: neither source dominates, so both are measured well in one run.
Distance dependence:

| distance from face | Ω/4π | interacting | ratio vs La-138 |
|---|---|---|---|
| 100 mm | 3.28% | 974 /s | 1.79 |
| 150 mm | 1.54% | 457 /s | 0.84 |
| **200 mm** | **0.883%** | **263 /s** | **0.48** |
| 300 mm | 0.398% | 118 /s | 0.22 |
| 500 mm | 0.145% | 43 /s | 0.08 |

**A distance scan is the strongest test available here.** The source rate must fall as 1/d²
while the La-138 rate stays flat — it is internal and cannot depend on where the source sits.
If the intrinsic rate moves with distance, something is wrong with the source classification
or the geometry, not with the activity calculation. Run three or four distances and check
both trends.

## Statistics

Livetime = `maxEvents` × 1695 ns; 100,000 events = 0.1695 s. At `activityScale : 100` that
gives ~4,500 source pulses and ~9,200 La-138 pulses — enough for a few-percent rate
comparison. Raise both scales together for more.

## Source classification

`LaBrTruthPulses` separates the two by the **start z of the traced primary**, threshold
`photonGunMaxZ : 40850.0`. The external source sits upstream of the crystal front face at
40853.4 mm; La-138 vertices are inside the crystal. The 3.4 mm margin is comfortable for any
realistic source position, but a source placed closer than that to the face would be
misclassified as intrinsic — raise the threshold if you ever do that.

Note the reused label: the external source appears as `source == 0`, the `photonGun` enum
value in `STMMCPulse`. No data-product change was needed for that.

## Caveats

- **This is the Mu2e geometry, not a bench.** The simulated source sits inside the STM
  enclosure with the SSC collimator and shielding present. A real bench measurement has none
  of that. The absolute efficiency will therefore differ; the **ratio** and the **1/d² trend**
  are the robust comparisons, not the raw efficiency.
- **Attenuation between source and crystal is whatever the Mu2e geometry contains.** Check
  the digitizer's `endJob` loss breakdown to see how many steps are actually reaching the
  crystal.
- Beta particles from the source are not emitted — sealed sources absorb them in their
  encapsulation. Correct for Cs-137, Co-60 and Mn-54; for a bare source it would not be.
- The intrinsic peak efficiency used in the prediction above (~0.55) is a literature estimate
  for a 3″×3″ LaBr₃, not a measurement. The simulation computes the real one; the analytic
  number is only there as an order-of-magnitude cross-check.
