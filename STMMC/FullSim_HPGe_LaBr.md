# Full STM simulation: HPGe + LaBr together (photon gun + intrinsic La-138)

A single simulation that exercises **both** STM detectors so their energy spectra can be
plotted together. Photon-gun signal goes into each detector; intrinsic La-138 self-
scintillation is added in the LaBr crystal. Read this together with `La138_STATUS.md`
(the La-138 / Method C entry point).

## What it does

`fcl/FullSimToHPGeAndLaBrWaveforms.fcl` wires:

- **genHPGe** (`PhotonGun`) - 1.809 MeV gamma fired along +z from the HPGe aperture.
- **genLaBr** (`PhotonGun`) - 1.809 MeV gamma fired along +z from the LaBr aperture.
- **la138** (`La138Source`, Method C) - intrinsic La-138 decay products emitted from uniform
  vertices inside the LaBr crystal. La-138 is intrinsic to LaBr3 only (germanium has no
  lanthanum), so there is deliberately no HPGe intrinsic source.
- **genCombined** (`GenParticleMerger`) - concatenates the three generator collections.
- **g4run** - Geant4 EM transport, `turnOnRadioactiveDecay: false` (Method C emits the
  products directly), reading `genCombined`.
- **DigiHPGe** (`HPGeWaveformsFromStepPointMCs`) and **DigiLaBr**
  (`LaBrWaveformsFromStepPointMCsTest`) - both read the shared `g4run:STMDet` collection and
  self-select their own crystal by the x cut (HPGe keeps x < -3904, LaBr keeps x > -3925),
  so one collection feeds both digitizers.
- **TruthLaBr** (`LaBrTruthSpectrum`) - Geant4 truth spectrum for the LaBr crystal, split
  into photon-gun vs La-138 (LaBr-only; see "Open" below).

Both digitizers write a per-ADC-bin TTree (`DigiHPGe/ttree`, `DigiLaBr/ttree`) with branch
`chargeCollected`. Summing `|chargeCollected|` over an event gives that detector's deposited
charge; the plot macro converts it to keV.

## HPGe module promotion (one-time integration step, done)

The improved HPGe digitizer lived in `src/HPGeWaveformsFromStepPointMCs_module_test.cc`,
which **was not built** (CMake builds only `*_module.cc`, and an art plugin file must end in
`_module.cc`). The old `_module.cc` had the wrong crystal centre (z off by ~87 mm), which
rejected ~99% of in-crystal HPGe steps. The improved file was therefore **copied over**
`HPGeWaveformsFromStepPointMCs_module.cc` so the fixed-geometry / Fano-factor version is the
one built. `_module_test.cc` is now redundant (kept as the reference copy, but it is now
DIVERGENT from `_module.cc` - see the RNG change below; treat `_module.cc` as canonical).

**RNG change (required to make it link).** The improved module originally seeded its noise
via `SeedService` (`createEngine(ServiceHandle<SeedService>()->getSeed())` + `CLHEP::RandGaussQ`).
On the gpvm this would **not link**: `undefined reference to mu2e::SeedService::getSeed()`, and
adding `Offline::SeedService` to `CMakeLists.txt` did NOT fix it (the long-standing STMMC
SeedService linking quirk - the LaBr module hit the same wall). It was therefore converted to
`std::default_random_engine _randomGen` + `std::normal_distribution<double> _noiseGauss{0,1}`
(called as `_noiseGauss(_randomGen)`), mirroring the LaBr digitizer. The includes
`RandGaussQ.h`, `RandomNumberGenerator.h`, `SeedService.hh` were removed and the `CMakeLists.txt`
`Offline::SeedService` addition was reverted (net no CMake change). **Do not re-introduce
SeedService into any STMMC module - use stdlib RNG.**

> Side effect: the improved module does `produces<STMWaveformDigiCollection>("HPGe")`
> (instance name `"HPGe"`), whereas the old one produced the empty instance. The prolog reco
> tags that read the raw digi (`HPGeDigitization.concatenation.STMWaveformDigisTag` and
> `STMMCAnalysis.ZS.HPGe...uSpill`) were updated `DigiHPGe:` -> `DigiHPGe:HPGe`, and the LaBr
> equivalent `STMMCAnalysis.ZS.LaBr.stmWaveformDigisTag` `DigiLaBr:` -> `DigiLaBr:LaBr` (the
> LaBr digitizer likewise produces the `"LaBr"` instance). `HPGeReco.fcl` now also sets
> `FanoFactor : @local::HPGeDigitization.FanoFactor`, and `HPGeDigitization.FanoFactor : 0.1`
> was added to `prolog.fcl` (it was referenced but undefined).

## Build / run / plot

```bash
# Build is required once because the HPGe C++ module changed (promotion above).
muse setup
muse build -j8

mu2e -c Offline/STMMC/fcl/FullSimToHPGeAndLaBrWaveforms.fcl -n 1000
root -l -b -q 'Offline/STMMC/src/plotHPGeAndLaBrSpectra.C("nts.FullSimHPGeLaBr.root")'
```

Outputs:
- `sim.FullSimHPGeLaBr.art` - art event data
- `nts.FullSimHPGeLaBr.root` - `DigiHPGe/ttree`, `DigiLaBr/ttree`, `TruthLaBr/ttree`
- `HPGeLaBr_EnergySpectra.pdf` - both detectors overlaid on one canvas
- `HPGeLaBr_EnergySpectra_panels.pdf` - the two spectra side by side
- (optional) `root -l -b -q 'Offline/STMMC/src/plotLaBrWaveform.C("nts.FullSimHPGeLaBr.root")'`
  and `plotLa138TruthSpectrum.C` still work on the LaBr/truth trees in the same file.

## Build / bring-up troubleshooting (lessons learned on the gpvm)

The user edits locally and copies files to the gpvm (`/exp/mu2e/app/users/sclarke/Simulation`)
to build/run; **nothing builds locally**. The bring-up of this work hit, in order:

1. **Runtime "module misconfigured" (missing `StepPointMCsTag`/`FanoFactor`, requires
   `StepPointMCsTagEle/Mu/1809`).** This is the OLD module's config signature - it means the
   OLD plugin `.so` was loaded. Root cause was a **failed link** (next item): a build that
   errors leaves the previous `.so` in place, so the runtime keeps loading old code even though
   the source is new. **Always confirm the library actually LINKED, not just that the build ran
   and the source looks right.** Quick checks on the gpvm:
   - `grep -n "FanoFactor\|crystalCentreX =" Offline/STMMC/src/HPGeWaveformsFromStepPointMCs_module.cc`
     (expect 2x FanoFactor and `-3986.30`)
   - `touch <the _module.cc>; muse build -j8 2>&1 | tee /tmp/bld.log` then
     `grep -iE "HPGeWaveforms|error|FAILED" /tmp/bld.log` - watch the `.so` LINK line complete.
2. **`undefined reference to mu2e::SeedService::getSeed()`** at link - see the RNG change above.
   Fixed by switching to stdlib RNG, not by CMake.

If a future change to an STMMC module pulls in a new dependency, expect link surprises; match
the library list of an already-working STMMC module first.

## Energy estimator (how the spectra are built)

`chargeCollected` is filled per ADC bin as the per-bin charge increment, so its sum over an
event equals the total charge that detector collected (one deposit's full charge, modulo
window-edge carry-over). Per event:

- **HPGe**: charge is in Ge electron-hole pairs. `E[keV] = |sum| * 2.96 / 1000`
  (epsilonGe = 2.96 eV/pair). Carries the Fano fluctuation if `FanoFactor > 0`.
- **LaBr**: charge is in photoelectrons. `E[keV] = |sum| / (63000*0.35) * 1000`.
  Carries Poisson photostatistics + PMT excess-noise smearing.

This is the pre-ADC charge estimator (same as `plotLaBrWaveform.C`): it includes detector
resolution but not ADC quantization or electronic noise.

## Tuning

- `la138.activityScale` (default 1000) multiplies the physical ~544 Bq (~9e-4 decays/event).
  At 1000 -> ~0.9 decays/event in LaBr (visible intrinsic peaks). Set to 1 for physical
  normalization. La-138 gammas can Compton-scatter into the HPGe (real cross-talk); a high
  activityScale inflates that, so lower it for a clean HPGe signal study.
- `maxEvents` / `-n` sets statistics; Geant4 is the slow step.

## Current status (2026-06-18)

- HPGe module promoted + converted to stdlib RNG; **expected to link cleanly now**. Confirm the
  next `muse build -j8` finishes the `libmu2e_STMMC_HPGeWaveformsFromStepPointMCs_module.so`
  LINK line with no `undefined reference`.
- prolog `FanoFactor`, reco instance-tag fixes, and `HPGeReco.fcl` Fano line are all in.
- **Plotting not yet confirmed working.** `plotHPGeAndLaBrSpectra.C` was run on the gpvm and
  produced NO console output and no PDFs. The macro prints `Building spectra from ...` as its
  first line (and `Cannot open ...` on a bad file), so zero output means the macro that ran was
  NOT this one - almost certainly the macro file was not copied to the gpvm (or a stale/empty
  copy is there). Next session, verify on the gpvm before anything else:
  - `wc -l Offline/STMMC/src/plotHPGeAndLaBrSpectra.C` (expect ~165) and
    `grep -c "Building spectra from" Offline/STMMC/src/plotHPGeAndLaBrSpectra.C` (expect 1).
    If 0/empty -> recopy the macro (no build needed, it is interpreted ROOT).
  - `root -l -b -q -e 'auto f=TFile::Open("nts.FullSimHPGeLaBr.root"); if(f) f->ls();'` -
    expect `DigiHPGe`, `DigiLaBr`, `TruthLaBr` directories (each with `ttree`). If missing, the
    digitizers did not write TTrees (check the run finished and `makeTTree:true`).

## Open / next steps

1. **HPGe truth spectrum**: `LaBrTruthSpectrum` is LaBr-specific (its cylinder cut assumes the
   LaBr crystal's identity rotation). The HPGe crystal is rotated (~45 deg about y) and has a
   bore hole, so an HPGe truth analyzer must replicate the digitizer's rotation + bore
   handling. Not done here. The HPGe *digitized* spectrum is produced and plotted.
2. **End-to-end gpvm validation**: confirm both photopeaks land at 1809 keV and the LaBr
   intrinsic peaks (1468/1436/789 keV) appear, with sane counts.
3. The `_module_test.cc` reference copy can be deleted once the promotion is confirmed good
   (it is now divergent from `_module.cc`; `_module.cc` is canonical).

## Files this work touches (copy to gpvm, then `muse build -j8` since C++ changed)

- `Offline/STMMC/src/HPGeWaveformsFromStepPointMCs_module.cc` (promoted geometry + stdlib RNG)
- `Offline/STMMC/fcl/FullSimToHPGeAndLaBrWaveforms.fcl` (new combined pipeline)
- `Offline/STMMC/src/plotHPGeAndLaBrSpectra.C` (new overlay plot macro; ROOT-interpreted)
- `Offline/STMMC/fcl/prolog.fcl` (added `HPGeDigitization.FanoFactor : 0.1`; reco tags ->
  `DigiHPGe:HPGe` / `DigiLaBr:LaBr`)
- `Offline/STMMC/fcl/HPGeReco.fcl` (added `FanoFactor` line)
- `Offline/STMMC/CMakeLists.txt` is net-unchanged (SeedService add was reverted)
