# Comparing the modular LaBr pipeline with the old job files

How to check that `LaBrPipeline*.fcl` (with `STMParticleSource`) reproduces the old one-off
LaBr jobs.

- **Gun test: expect an exact, bin-for-bin match** (ratio = 1 everywhere). Neither gun draws a
  random number (fixed point, direction, energy, t = 0), g4run gets the same SeedService seed
  in both jobs, and the digitizer's `std::default_random_engine` is never seeded, so it
  replays the same sequence in every job. Same primaries + same seeds = the same events.
- **Poisson + La-138 test: expect statistical agreement only.** Here the generators do draw
  random numbers (Poisson counts, times, decay vertices, branches, beta energies), with
  different code (PoissonPhotonGun + La138Source vs STMParticleSource + LaBrIntrinsicSource),
  so the primaries differ and every event after that differs. The test is whether the
  per-event spectra and rates agree within statistical errors.

Two paired tests cover all the new pieces that have an old counterpart:

| Test | Old job | New job | Checks |
|---|---|---|---|
| Gun | `PhotonGunToLaBrWaveform.fcl` | `LaBrPipeline.fcl` (default) | gun mode, `single` emission, pipeline wiring |
| PoissonLa138 | `PoissonPhotonGunPlusLa138ToLaBrWaveform.fcl` | `NewPoissonLa138.fcl` (below) | `poisson` emission, `intrinsic` source, merging, truth |

In each pair the detector physics is the same: the same digitizer module and parameters,
`QGSP_BERT_EMZ`, `STMDet` only, the same beam point (LaBr aperture axis, z = 40404 mm, +z),
1.809 MeV photons, a 1695 ns window and 10000 events. The only difference is the generator.

---

## 0. Set up and build

On the gpvm, in `/exp/mu2e/app/users/sclarke/Simulation`:

```bash
muse setup
muse build -j8 2>&1 | tail -5
find build -name "*STMParticleSource*.so" -o -name "*STMParticleDump*.so"
```

Both `.so` files must be listed. A build that fails to link leaves the old plugins in place, so
a missing file means the new modules were not built. Files needed on the gpvm:
`src/STMParticleSource_module.cc`, `src/STMParticleDump_module.cc`, `src/compareLaBrSpectra.C`,
`fcl/LaBrPipelineProlog.fcl`, `fcl/LaBrPipelineJob.fcl`, `fcl/LaBrPipeline.fcl`, `CMakeLists.txt`.

Each job runs in its own directory. Several jobs write output files with the same names, and
`plotLaBrWaveform.C` always writes the same PDF names.

```bash
mkdir -p cmp/old_gun cmp/new_gun cmp/old_poisson cmp/new_poisson
```

---

## 1. Gun test

```bash
(cd cmp/old_gun && mu2e -c Offline/STMMC/fcl/PhotonGunToLaBrWaveform.fcl -n 10000)
(cd cmp/new_gun && mu2e -c Offline/STMMC/fcl/LaBrPipeline.fcl)
```

The old file sets `maxEvents : 100`, so `-n 10000` is required. The new default is
`LaBrRun.nEvents : 10000`.

```bash
root -l -b -q 'Offline/STMMC/src/compareLaBrSpectra.C+("cmp/old_gun/nts.PhotonGunToLaBrWaveform.root","cmp/new_gun/nts.LaBrPipeline.root","Gun")'
```

The old gun job has no `TruthLaBr` analyzer, so only the digitized spectrum is compared.
The macro prints "truth comparison skipped" for this test. That is expected.

---

## 2. Poisson + La-138 test

The old job fires Poisson(0.0339) photons per event (20 kHz) with `PoissonPhotonGun`, plus
La-138 decay products at `activityScale : 1000` from `La138Source`. To match it, the new job uses the
`poisson` emission (20 kHz) and `LaBrIntrinsicSource` with the Ac-227 chain off, at the same
scale. Both La-138 generators use the same decay data (`inc/LaBrIntrinsicDecays.hh`) and the
same activity computed from the crystal composition.

```bash
(cd cmp/old_poisson && mu2e -c Offline/STMMC/fcl/PoissonPhotonGunPlusLa138ToLaBrWaveform.fcl)

cat > cmp/new_poisson/NewPoissonLa138.fcl <<'EOF'
#include "Offline/STMMC/fcl/LaBrPipelineProlog.fcl"
BEGIN_PROLOG
LaBrRun.processName : "NewPoissonLa138"
LaBrRun.nEvents : 10000
LaBrRun.primary.emission : @local::LaBrPipeline.Emission.poisson
LaBrRun.sources : [ primary, intrinsic ]
LaBrRun.intrinsic.includeAc227 : false
LaBrRun.intrinsic.activityScale : 1000.0
LaBrRun.histFile : "nts.NewPoissonLa138.root"
LaBrRun.artFile : "sim.NewPoissonLa138.art"
END_PROLOG
#include "Offline/STMMC/fcl/LaBrPipelineJob.fcl"
EOF
(cd cmp/new_poisson && mu2e -c $PWD/NewPoissonLa138.fcl)
```

```bash
root -l -b -q 'Offline/STMMC/src/compareLaBrSpectra.C+("cmp/old_poisson/nts.PoissonPhotonGunPlusLa138.root","cmp/new_poisson/nts.NewPoissonLa138.root","PoissonLa138")'
```

---

## 3. Optional: waveform shapes

```bash
for d in old_gun new_gun old_poisson new_poisson; do
  (cd cmp/$d && root -l -b -q "../../Offline/STMMC/src/plotLaBrWaveform.C(\"$(ls nts.*.root)\")")
done
```

Compare `cmp/*/LaBr_IndividualPulses.pdf` by eye. The pulse shape comes from the digitizer,
so old and new should look the same. A difference points to the digitizer configuration, not
the generator.

---

## 4. What the macro produces

`compareLaBrSpectra.C(old, new, tag)` normalises both runs per event (micropulse) and writes:

- `LaBrCompare_<tag>_Digitized.pdf`: the digitized energy per event, as in `plotLaBrWaveform.C`:
  the sum of |chargeCollected| times 1000/(63000 x 0.35) keV, in 250 bins from 0 to 2000 keV.
  Old is solid blue, new is dashed red, with a new/old ratio panel below.
- `LaBrCompare_<tag>_Truth_External.pdf` and `..._Truth_La138.pdf`: the Geant4 truth energy
  per decay from `TruthLaBr/decays`, split by source (0 = external, 1 = La-138).

It prints:

```
=== <tag>: <old> (old) vs <new> (new)
  events: old N, new N
  events with a deposit: old x, new y per event
  mean digitized energy of hit events [keV]: old ..., new ...
  digitized: entries ... | chi2 prob ... | KS prob ...
  truth decays per event (+- stat), by source:
    external: old ... +- ..., new ... +- ...
    La-138:   old ... +- ..., new ... +- ...
```

---

## 5. Pass criteria

| Quantity | Expected | Worth investigating |
|---|---|---|
| events | 10000 in both | a different count: the event counts or `-n` don't match |
| events with a deposit / event | equal within ~3 sigma, sigma ~ sqrt(n)/N | a systematic offset larger than that |
| mean digitized energy | equal within a few keV | a shift of tens of keV |
| peak positions | 1809 photopeak, ~1585 Compton edge, 1298 single escape, 787 double escape (gun); plus 1436/1468 and the 789-1052 continuum (La-138) | any peak moved or missing |
| ratio panel | scatter around 1, wider where counts are low | a trend or step with energy |
| chi2 prob, KS prob | roughly uniform in 0-1 over repeated seeds | consistently < 0.01 |
| La-138 decays / event | the same in old and new; about 505 Bq x 1000 x 1695e-9 = 0.86, slightly less because only decays that deposit energy are counted | disagreement beyond the printed +- errors |
| external decays / event | the same in old and new (0.0339 thrown per event, times the fraction that reach the crystal) | disagreement beyond the printed +- errors |

Notes on reading these numbers:
- Applying the KS test to binned histograms gives probabilities that are biased high, so treat
  chi2 as the main test and KS as a sanity check.
- One pass at p ~ 0.05 means nothing. To be sure, rerun the new job with a different seed and
  check that the p-values move around. Set `LaBrRun.baseSeed : 9` in the BEGIN_PROLOG block, or
  add `services.SeedService.baseSeed : 9` at the end of the old fcl.
- The old job `PoissonPhotonGunPlusLa138` counts decays with the older `La138Source`. If a
  La-138 difference appears only in the "external"/"La-138" split and not in the total, check
  source classification (GenId) rather than the physics.

## 6. If a test fails

1. Check that the new job actually used the new modules. The log should contain an
   `STMParticleSource` summary line (`emitted ... particles in ... events`).
2. Gun test: `emitted` should equal 10000 (one photon per event).
3. Poisson test: `emitted` / events should be 0.0339 +- 0.0018 for 10000 events.
4. Compare the job configurations directly (`--debug-config` writes the fully expanded configuration and exits):
   ```bash
   (cd cmp/new_poisson && mu2e -c $PWD/NewPoissonLa138.fcl --debug-config new.cfg)
   (cd cmp/old_poisson && mu2e -c Offline/STMMC/fcl/PoissonPhotonGunPlusLa138ToLaBrWaveform.fcl --debug-config old.cfg)
   diff <(grep -A30 "DigiLaBr" cmp/old_poisson/old.cfg) <(grep -A30 "DigiLaBr" cmp/new_poisson/new.cfg)
   ```
   The only expected digitizer difference is `alphaQuenchingFactor` / `alphaQuenchingSlope`,
   which only affect alphas and so are not tested here.
