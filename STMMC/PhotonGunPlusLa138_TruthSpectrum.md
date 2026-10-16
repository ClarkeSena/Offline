# PhotonGunPlusLa138 truth-spectrum work

Status note for a future Claude session. Covers the `PhotonGunPlusLa138` combined
simulation and the per-source Geant4 **truth** energy spectrum built on top of it.

## Goal

`PhotonGunPlusLa138ToLaBrWaveform.fcl` runs two generators into one g4run:
- **PhotonGun** — the 1.809 MeV signal gamma, fired straight along +z (`px=py=0`) from the
  STM aperture at (x,y,z) = (LaBr.x, LaBr.y, 40404), upstream of the crystal.
- **La138Source** — intrinsic La-138 decay products (Method C) emitted from uniform
  vertices *inside* the LaBr crystal (z ≈ 40891 ± 38 mm).

Both deposit into the same `g4run:STMDet` StepPointMC collection, so the digitized waveform
mixes them. The user wanted the **truth** (Geant4-deposited) energy spectrum split into
three plots: combined, photon-gun-only, La138-only.

## What was added

### 1. `src/LaBrTruthSpectrum_module.cc` (new EDAnalyzer)
Reads `g4run:STMDet`, sums per-event ionizing energy into three TTree branches
(`edepTotal`, `edepGun`, `edepLa138`, all in **keV**) in `TruthLaBr/ttree`.

**Source classification (the crux):** each step is traced to its true primary, classified
by the primary's **start z** vs `photonGunMaxZ` (fhicl, default 40700, sitting in the clean
gap between the aperture at 40404 and the crystal at ~40891). Below → photon gun, above →
La-138.

> **Important gotcha:** `SimParticle::originParticle()` does **NOT** walk to the primary.
> It only recurses for the special "self-parent" case
> (`selfParent() = hasParent() && creationCode == mu2ePrimary`). For an ordinary secondary
> (a Compton/photo electron, creationCode `compt`/`phot`) it returns *the secondary itself*,
> whose start vertex is inside the crystal — so everything got misclassified as La-138.
> The module instead walks `parent()` explicitly up to the no-parent primary. Do not
> "simplify" this back to `originParticle()`.

**Selection cuts** (match the digitizer so truth ↔ digitized are comparable):
- LaBr-side x cut: `x > -3925` (rejects the HPGe crystal, which shares the collection).
- LaBr crystal-cylinder bounds (centre (-3863.36, 0, 40891.5), L=76.2, R=38.1, tol 0.1).
- No time cut.
- Skips `ionizingEdep() <= 0` steps (gammas/neutrinos — zero energy, no effect).

(History: cuts were briefly removed then restored at the user's request. Current state =
x cut + cylinder bounds, no time cut.)

CMake: `cet_build_plugin(LaBrTruthSpectrum ...)` linking `MCDataProducts` +
`art_root_io::TFileService_service`.

### 2. `fcl/PhotonGunPlusLa138ToLaBrWaveform.fcl`
- Added `analyzers.TruthLaBr` (module `LaBrTruthSpectrum`, `StepPointMCsTag: "g4run:STMDet"`,
  `photonGunMaxZ: 40700`) and a `truth_path` in `end_paths`.
- **`la138.activityScale` was uncommented and set to 1000.0.** It had been left commented
  (`# activityScale : 1.0`), so editing that line did nothing — the module fell back to the
  default 1.0 (physical ~540 Bq → 9e-4 decays/event). At 1000 → ~0.9 decays/event (≈1 per
  event, keeps the per-event sum clean). `la138` has `verbosityLevel: 1`; the startup log
  prints `activity = ... Bq (scaled)` and `mean decays per ... event` — check it to confirm
  the scale took.

### 3. Plotting macros (in `src/`)
- `plotLa138TruthSpectrum.C` — reads `TruthLaBr/ttree`, makes 3 PDFs:
  `LaBr_TruthSpectrum_{Combined,PhotonGun,La138}.pdf`. Fixed binning 250 bins, 0–2000 keV.
- `plotLaBrWaveform.C` — existing macro, reads `DigiLaBr/ttree` (full waveform, individual
  pulses, **digitized** energy spectrum). Its energy spectrum was changed from auto-range to
  the **same fixed 250 bins / 0–2000 keV** so it overlays the truth spectrum bin-for-bin.

## Build / run / plot

```bash
muse build -j8          # needed only after C++ (module) changes, not for fcl edits
mu2e -c Offline/STMMC/fcl/PhotonGunPlusLa138ToLaBrWaveform.fcl
root -l -b -q 'Offline/STMMC/src/plotLa138TruthSpectrum.C("nts.PhotonGunPlusLa138.root")'
root -l -b -q 'Offline/STMMC/src/plotLaBrWaveform.C("nts.PhotonGunPlusLa138.root")'
```

Both trees live in the same `nts.PhotonGunPlusLa138.root` (one TFileService output).

## Physics interpretation (settled during this session)

**Combined truth spectrum is correct, including the continuum below 1809.** A 1809 keV gamma
in a finite crystal gives, besides the photopeak: Compton continuum + edge at ~1585 keV
(2E²/(m_ec²+2E)), single-escape at 1298 keV (1809−511), double-escape at 787 keV
(1809−1022, under the La 789 line), plus partial deposits from scattering in passive
material between aperture and crystal. None of this is a bug.

**Truth vs digitized differ by detector resolution, not calibration.** Same pe/MeV constants
(`63000 × 0.35 = 22050 pe/MeV`) → same peak positions and integrals. The digitized photopeak
is a Gaussian; the truth one is a delta. Smearing is introduced in the digitizer
`LaBrWaveformsFromStepPointMCsTest::depositCharge()`:
1. **Poisson photostatistics** — `N_pe = Poisson(Edep × 22050)`. σ_E/E = 1/√N̄ (shrinks as
   1/√E). Counting statistics of photoelectron production.
2. **PMT excess noise factor F = 1.27** — `signal_pe += Gauss(0, √((F−1)·N_pe))`. Dynode-chain
   gain jitter; σ_E/E = √(F/N̄). Tuned to reproduce the measured 2.2% FWHM at 662 keV.
   At 1809 keV this gives ~1.3% FWHM ≈ 24 keV.

Other (minor) truth↔digitized differences: ADC quantization; inter-micropulse carry-over
(`microspillBufferLengthCount = 2`) can move a pulse tail into the next event's window.
`NoiseSD` is 0.0 in this fcl, so electronic noise contributes nothing.

## Known caveats / possible next steps
- The module fills **per event** (summed). At `activityScale` high enough for >1 La-138
  decay/event, multiple decays pile up within an event and smear the 1436/789 keV lines.
  For more stats with clean lines, raise `maxEvents` at `activityScale ≈ 1000` (≤1
  decay/event), or refactor to **per-primary** accumulation (group each step's energy by the
  primary it traces to, fill once per primary). The per-primary refactor was offered but not
  done.
- La-138 decay-scheme constants in `La138Source_module.cc` are marked
  "VERIFY/REFINE before publication" (EC fraction, energies, X-ray probs, beta shape with
  Fermi function omitted).
- The crystal geometry z (40891.5) carries the unexplained −101.6 mm empirical offset noted
  in CLAUDE.md (`STMGeometryDump` reports z 101.6 mm high). Root cause still unknown.
