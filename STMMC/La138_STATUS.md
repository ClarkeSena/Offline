# La-138 intrinsic activity + truth spectrum — STATUS / session handoff

**Last updated:** 2026-06-17. **Read this first**, then the two detailed docs it points to.

This is the single entry point for the LaBr3 intrinsic-radioactivity (¹³⁸La self-
scintillation) work and the truth-spectrum tooling built around it. It records what is
done, what is verified, what is still open, and the gotchas not to re-discover.

---

## TL;DR state

- **Method C is implemented and standalone-validated.** `La138Source` emits the ¹³⁸La
  decay *products* directly (no ion, no G4 radioactive decay), at the Poisson rate set by
  the crystal's lanthanum content. G4 transports them, so the intrinsic spectrum features
  emerge from transport.
- **Truth-spectrum tooling is implemented** (`LaBrTruthSpectrum` analyzer +
  `plotLa138TruthSpectrum.C`): builds the Geant4 truth energy spectrum split into
  photon-gun signal vs La-138 intrinsic, for bin-for-bin comparison with the digitized
  spectrum.
- **What has NOT been independently confirmed by me:** that a full gpvm build + run
  produces the expected peak counts (the physics/logic is validated standalone; the
  end-to-end Geant4 run is the remaining check). See "Open / next" below.

---

## File map (all under `Offline/STMMC/`)

**Modules (`src/`)**
- `La138Source_module.cc` — Method C generator: per Poisson-drawn decay, emits products at
  a uniform crystal vertex sharing position+time. EC 66.4% -> 1.436 MeV γ + (71%) Ba K
  X-ray (0.0322/0.0364 MeV); β⁻ 33.6% -> 0.789 MeV γ + electron from the allowed β
  spectrum (endpoint 0.263 MeV). Activity computed from geometry/material (~544 Bq).
- `GenParticleMerger_module.cc` — concatenates GenParticleCollections so PhotonGun +
  La138Source feed one g4run (point `g4run.inputs.primaryTag` at the merged collection).
- `LaBrTruthSpectrum_module.cc` — EDAnalyzer; traces each `g4run:STMDet` step to its true
  PRIMARY (walking `parent()` manually — see gotcha) and bins per-event deposited energy
  into `TruthLaBr/ttree` branches `edepTotal`/`edepGun`/`edepLa138` (keV). Split by primary
  start-z vs `photonGunMaxZ` (default 40700).
- `LaBrWaveformsFromStepPointMCsTest_module.cc` — the LaBr digitizer (fitted pulse shape,
  peak-height calibration `peakKeVPerADC`).
- `BinaryWaveformWriter_module.cc` — raw int16 waveform stream output.
- `PoissonPhotonGun_module.cc` — Poisson-in-time photon gun (20 kHz signal stream).

**Plot macros (`src/`)**
- `plotLaBrWaveform.C` — waveforms + per-event energy spectrum (250 bins 0-2000 keV).
- `plotLa138TruthSpectrum.C` — 3 PDFs (Combined / PhotonGun / La138) from `TruthLaBr/ttree`,
  same binning for overlay with the digitized spectrum.

**FHiCL (`fcl/`)**
- `La138OnlyToLaBrWaveform.fcl` — La-138 ALONE (no gun). Use to see the intrinsic spectrum
  cleanly. `activityScale : 100` default.
- `PhotonGunPlusLa138ToLaBrWaveform.fcl` — combined signal + intrinsic (plain PhotonGun,
  1 photon/event at t=0; La138Source Poisson). `turnOnRadioactiveDecay: false`.
- `PoissonPhotonGunToLaBrWaveform.fcl` — 20 kHz Poisson signal stream.
- `PhotonGunToLaBrWaveform.fcl` — single-photon signal test.

**Docs (`Offline/STMMC/`)**
- `La138_MethodC_implementation.md` — how Method C works + was built (full physics + code).
- `PhotonGunPlusLa138_TruthSpectrum.md` — the truth-spectrum analyzer + interpretation.
- `FullSim_HPGe_LaBr.md` — combined HPGe+LaBr sim (photon gun into both + La-138 in LaBr),
  plotting both energy spectra together, the HPGe module promotion + stdlib-RNG fix, and the
  gpvm build/plot troubleshooting. Read it for any two-detector / HPGe-integration work.
- `La138_STATUS.md` — THIS file.

---

## How to run

```bash
muse build -j8

# Intrinsic spectrum alone (cleanest validation; raise activityScale, many events)
mu2e -c Offline/STMMC/fcl/La138OnlyToLaBrWaveform.fcl -n 20000
root -l -b -q 'Offline/STMMC/src/plotLaBrWaveform.C("nts.La138Only.root")'

# Combined signal + intrinsic, with truth spectrum by source
mu2e -c Offline/STMMC/fcl/PhotonGunPlusLa138ToLaBrWaveform.fcl -n 14000
root -l -b -q 'Offline/STMMC/src/plotLa138TruthSpectrum.C("nts.PhotonGunPlusLa138.root")'
```

`activityScale` (in the `la138` fcl block) multiplies the physical ~544 Bq. Physical rate
is ~9.2e-4 decays/event, so for visible statistics use 100 (clean, low pile-up, run many
events) or ~1000 (faster, ~0.9/event but ~22% multi-decay pile-up). Set to 1 for physical
normalization in a real combined run.

## Expected counts (predict vs observe)

For `La138Only`, `activityScale=100`, `N` events: total decays `D = N * 100 * 544 * 1695e-9
= N * 0.0922`. Ideal (full-absorption) split: 1468 keV sum peak `0.471*D`, 1436 keV escape
peak `0.193*D`, 789-1052 keV β continuum `0.336*D`. Measured peaks are LOWER by the γ
photopeak efficiency (Compton escape, ~0.4-0.6) — that ratio measures the efficiency. A
feature at ZERO is the red flag that products are not landing in the crystal.

---

## Open / next steps

1. **End-to-end gpvm validation** (the main remaining check): confirm `STMDet`
   StepPointMCs appear with the right particles/energies and vertices inside the crystal,
   and that peak counts roughly match the predictions above.
2. **Auger simplification** — for the 29% of EC with no fluorescence X-ray the code emits
   nothing for the ~32 keV K-shell energy; physically those are Auger electrons (absorbed
   locally), so they should land at 1468 not 1436. Current code over-populates 1436. Fix:
   always deposit the ~32 keV (as escapable X-ray with ω_K, else local electron).
3. **Per-primary accumulation** in `LaBrTruthSpectrum` (currently per-event): with >1 La
   decay per event the 1436/789 lines smear. Accumulate per primary instead.
4. **Fermi function** F(Z=58,E) on the β spectrum in `La138Source::sampleBetaKE()`.
5. **²²⁷Ac α bumps** — separate, bigger task: non-uniform contamination (radial profile,
   not flat sampling) + alpha quenching in the digitizer (Geant4 gives full α energy but
   the digitizer applies one lightYield to all particles, so α land too high).
6. For a realistic combined rate/timing run, swap the plain `PhotonGun` in the combined
   fcl for `PoissonPhotonGun` (meanPhotonsPerEvent 0.0339, maxTime 1695).

---

## Critical gotchas (do not re-discover)

- **Geant4 does not decay detector-material atoms** — radioactive decay acts only on ion
  TRACKS. That, plus the 10¹¹ yr half-life giving ~10⁻²⁵ in-window decay probability, is
  why injecting a La-138 ion + analogue `turnOnRadioactiveDecay` produces ZERO signal.
  Method C sidesteps this by emitting the products directly. Do NOT "fix" the combined fcl
  by turning radioactive decay back on.
- **`SimParticle::originParticle()` does NOT walk to the primary** — it only recurses for
  `selfParent()` (hasParent && creationCode==mu2ePrimary), so for ordinary secondaries it
  returns the secondary itself. `LaBrTruthSpectrum` walks `parent()` manually. Anything
  classifying steps by primary MUST do the same or it misclassifies everything.
- **Crystal z = 40891.50 mm**, NOT the 40993.1 from STMGeometryDump / `originInMu2e()`
  (101.6 mm = front tungsten too high, due to a reference-anchor mismatch between
  STMMaker.cc and constructSTM.cc). Hard-code 40891.50; do NOT query GeometryService.
- **GenParticle energies are in MeV** (1435.795 keV = 1.435795 MeV), position mm, time ns,
  Mu2e coordinates.
- **STMMC modules use `std::default_random_engine`**, not SeedService (linking limitation).
- **Internal ¹³⁸La gives NO sharp 789 keV peak** — the β almost always sums with the 789 γ,
  so it's a 789->1052 keV summed continuum. This is correct, not a bug, and distinguishes
  intrinsic ¹³⁸La from an external 789 keV source.
- **`activityScale` must be set, not just present** — `# activityScale` is a comment
  (no-op, defaults to 1.0). It must be uncommented to take effect.
- Module files must be named `*_module.cc`; FHiCL must be pure ASCII.

---

## One-paragraph summary

Intrinsic ¹³⁸La activity is simulated via Method C: `La138Source` emits the decay products
(γ + X-ray for EC, γ + β for β⁻) directly at the geometry-derived Poisson rate, G4
transports them, and the digitizer produces waveforms; `LaBrTruthSpectrum` +
`plotLa138TruthSpectrum.C` give the truth spectrum split by source for comparison with the
digitized spectrum. The decay-scheme logic is standalone-validated; the remaining work is
a full gpvm run to confirm counts, plus the refinements listed above (Auger, per-primary
binning, Fermi function, ²²⁷Ac + alpha quenching). See `La138_MethodC_implementation.md`
and `PhotonGunPlusLa138_TruthSpectrum.md` for the details.
