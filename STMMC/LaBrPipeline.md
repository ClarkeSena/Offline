# Modular LaBr pipeline

A single, configurable replacement for the one-off LaBr job files (`PhotonGunToLaBrWaveform`,
`PoissonPhotonGun*`, `*PlusLa138*`, ...). Those files are untouched and still work.

This is the quick reference. **`LaBrPipeline_Guide.md`** is the complete guide: step-by-step
internals, every parameter, recipes, and how to change the geometry.

```
[LaBrRun.sources] -> genCombined -> g4run -> DigiLaBr -> [TruthPulses, TagPileup] -> [TruthLaBr, BinaryWriter, fullOutput]
  primary     STMParticleSource   gun | spectrum | dataset, single | Poisson
  background  STMParticleSource   second slot, default = upstream dataset
  intrinsic   LaBrIntrinsicSource La-138 + Ac-227 chain self-activity
```

## Files

| File | Role |
|---|---|
| `src/STMParticleSource_module.cc` | New generator: one module, three orthogonal choices (mode / emission / beam) |
| `src/STMParticleDump_module.cc` | New analyzer: upstream VD StepPointMCs -> flat particle TTree (the "dataset") |
| `fcl/LaBrPipelineProlog.fcl` | Building blocks (`LaBrPipeline.*`) and the default run table (`LaBrRun`) |
| `fcl/LaBrPipelineJob.fcl` | Job body, built entirely from `LaBrRun`. Never edited |
| `fcl/LaBrPipeline.fcl` | **The file you edit / copy**: prolog include, `LaBrRun` overrides, job include |
| `fcl/LaBrPipeline_PoissonSpectrumPlusIntrinsic.fcl` | Example: line spectrum at 20 kHz + self-activity + .bin output |
| `fcl/LaBrPipeline_UpstreamDataset.fcl` | Example: 1809 signal + upstream dataset + self-activity |
| `fcl/DumpSTMParticles.fcl` | Makes the dataset tree from upstream art files |

## The three knobs (in the BEGIN_PROLOG block of `LaBrPipeline.fcl`)

```fcl
LaBrRun.primary          : @local::LaBrPipeline.Primary.gun        # gun | spectrum | dataset
LaBrRun.primary.emission : @local::LaBrPipeline.Emission.poisson   # single | poisson
LaBrRun.primary.emission.rateHz : 20000.0
LaBrRun.sources          : [ primary, intrinsic ]                  # drop intrinsic = no self-scintillation
```

Replace a whole table before changing keys inside it (FHiCL applies overrides in order).

### Input (`mode`)
- `gun`: `pdgId` (default 22) at kinetic energy `E` [MeV]. Any particle, e.g. `pdgId : 2112`.
- `spectrum`: `spectrum.type` =
  - `lines`: `lineEnergies` + relative `lineIntensities`
  - `histogram`: `binEdges` (N+1) + `binContents` (N, counts per bin), flat within a bin
  - `rootHist`: `fileName` + `histName` of a TH1, x axis in MeV
- `dataset`: resamples particles (with their own pdgId, position, momentum and time) from a TTree.
  `sampleUnit : "event"` keeps all particles of one upstream event together (correlated
  neutron + gamma, etc.); `"particle"` resamples single entries. `pdgIds`, `EMin`, `EMax` filter.
  `transform` maps the upstream plane onto the aperture (defaults: VD101 r=200 mm -> LaBr SSC
  aperture, transverse position and momentum scaled by 3.989/200, |p| kept), the same
  mapping as `ShiftVirtualDetectorStepPointMCs` in the production STMResampler, without its
  HPGe-absorber-shadow special case. Remove `transform` if the tree is already at the aperture.

Gun/spectrum start from `beam` (default: LaBr aperture axis at z=40404 mm, +z pencil beam).
`beam.radius` fills a disk, `beam.coneHalfAngle` opens a uniform cone (180 = isotropic).

### Emission (Poisson or not)
- `single`: exactly `nPerEvent` (1) per event at `time` 0: a classic particle gun.
- `poisson`: Poisson(`rateHz` x 1695 ns) per event (or give `nPerEvent` as the mean), times
  uniform in the micropulse. Dataset events keep their internal time structure from a random
  start; late particles spill into the next micropulse through the digitizer buffer.
- `timing : "dataset"` uses the upstream times folded into the window instead.

### Self-scintillation
`intrinsic` is `LaBrIntrinsicSource` (La-138 computed from composition, Ac-227 chain at
`LaBrIntrinsicActivity.ac227Chain`). `LaBrRun.intrinsic.includeAc227 : false` gives La-138 only;
`activityScale` boosts statistics (and pileup). Alpha quenching is applied in the digitizer and
both truth modules.

## Making a dataset from upstream simulation

```bash
mu2e -c Offline/STMMC/fcl/DumpSTMParticles.fcl -S upstream_VD101_files.txt   # -> nts.STMParticleDump.root
```
Tree `DumpParticles/particles`: eventId, run, subRun, event, virtualdetectorId, pdgId, trackId,
x, y, z [mm], px, py, pz [MeV/c], time [ns], E [MeV kinetic]. `DumpParticles/info` holds
`nEvents` (input events read). A tree from any other simulation works if it has `pdgId` (int),
`x y z px py pz time` (double) and, for `sampleUnit "event"`, `eventId` (Long64).

**Normalisation** (set `emission.nPerEvent`, it is not automatic):
- `sampleUnit "event"`: nPerEvent = upstream events per micropulse,
  e.g. (POT per micropulse) / (POT represented by one upstream event).
- `sampleUnit "particle"`: nPerEvent = that x (particles->GetEntries() / info.nEvents).
Particles with a PDG code missing from the particle table (e.g. recoil nuclei) are skipped and
counted in the log.

## Running

```bash
muse build -j8                                  # new modules
mu2e -c Offline/STMMC/fcl/LaBrPipeline.fcl      # -> sim/nts.LaBrPipeline.*
root -l -b -q 'Offline/STMMC/src/plotLaBrWaveform.C("nts.LaBrPipeline.root")'
```
All existing macros work on the output (same module labels: `DigiLaBr`, `TruthLaBr`, `TagPileup`).
Every STMParticleSource particle is `GenId::particleGun`, so the truth modules label it as the
external source (0), La-138 as 1 and Ac-227 as 2.

## Caveats
- Not yet built or run: compiled here only against stub headers.
- `TagPileup.nMerge` is tied to `LaBrRun.nEvents` (one batch). If you run with `-n`, set
  `LaBrRun.nEvents` to the same value.
- The digitizer micropulse is hard-coded at 1695 ns; `emission.windowTime` must stay 1695.
- The digitizer's random engine is never seeded (known issue, not yet fixed; see
  `STM_Simulation_Overview.md` §11). `LaBrRun.baseSeed` changes the generators and Geant4 but
  not the photostatistics / noise sequence, which is the same in every job.
