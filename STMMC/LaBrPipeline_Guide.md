# The modular LaBr pipeline: complete guide

This is the full reference for the modular LaBr3:Ce+Sr simulation. It covers how a job is
assembled, what happens to a particle from generation to the written ADC sample, every number
that links the stages, and how to change each part, including the geometry.
`LaBrPipeline.md` is the short version.

Contents
1. [Overview](#1-overview)
2. [Files and what each one does](#2-files-and-what-each-one-does)
3. [How a job is assembled (the FHiCL layer)](#3-how-a-job-is-assembled-the-fhicl-layer)
4. [One micropulse, step by step](#4-one-micropulse-step-by-step)
5. [The numbers that tie the stages together](#5-the-numbers-that-tie-the-stages-together)
6. [How to use it: recipes](#6-how-to-use-it-recipes)
7. [Changing the geometry](#7-changing-the-geometry)
8. [Outputs and the tools that read them](#8-outputs-and-the-tools-that-read-them)
9. [Parameters that must be kept in sync](#9-parameters-that-must-be-kept-in-sync)
10. [Validation status and known limitations](#10-validation-status-and-known-limitations)
11. [Troubleshooting](#11-troubleshooting)

---

## 1. Overview

```
 art event = one 1695 ns micropulse (EmptyEvent source, LaBrRun.nEvents events)

  +-----------+   +------------+   +-----------+
  | primary   |   | background |   | intrinsic |     generators named in LaBrRun.sources
  | STMParti- |   | STMParti-  |   | LaBrIntr- |     (each writes a GenParticleCollection)
  | cleSource |   | cleSource  |   | insicSrc  |
  +-----+-----+   +-----+------+   +-----+-----+
        +---------------+----------------+
                        v
                  genCombined  (GenParticleMerger: concatenates the collections)
                        v
                  g4run        (Mu2eG4: full Mu2e geometry, QGSP_BERT_EMZ, STMDet only)
                        |  -> g4run:STMDet StepPointMCs (HPGe + LaBr crystals) + SimParticles
            +-----------+-------------+------------------+
            v                         v                  v
        DigiLaBr                 TruthPulses        TruthLaBr (analyzer)
   (LaBrWaveformsFrom-         (one STMMCPulse      per-event + per-decay
    StepPointMCsTest)           per decay)           truth energies
            |                         v
            |                    TagPileup (merges all micropulses on one timeline,
            |                               tags overlapping pulses)
            v
   DigiLaBr:LaBr STMWaveformDigi (628 int16 ADC samples)
            +--> BinaryWriter  -> raw .bin stream
            +--> fullOutput    -> .art file
            +--> DigiLaBr/ttree (per-sample tree in the nts .root file)
```

The pipeline makes **three independent choices**, each set in one place:

| Choice | Setting | Options |
|---|---|---|
| What is thrown | `LaBrRun.primary` | `@local::LaBrPipeline.Primary.gun` / `.spectrum` / `.dataset` |
| How many and when | `LaBrRun.primary.emission` | `@local::LaBrPipeline.Emission.single` / `.poisson` |
| Crystal self-activity | `LaBrRun.sources` | include `intrinsic` or not |

Everything downstream of the generators (Geant4, digitizer, truth, outputs) is the same for
every choice.

---

## 2. Files and what each one does

All paths are under `Offline/STMMC/` unless shown otherwise.

### New in this pipeline
| File | Kind | Role |
|---|---|---|
| `src/STMParticleSource_module.cc` | EDProducer | The configurable generator (gun / spectrum / dataset x single / Poisson x beam shape) |
| `src/STMParticleDump_module.cc` | EDAnalyzer | Writes upstream virtual-detector crossings to the flat TTree that `dataset` mode reads |
| `fcl/LaBrPipelineProlog.fcl` | prolog | All building blocks (`LaBrPipeline.*`) and the default run table `LaBrRun` |
| `fcl/LaBrPipelineJob.fcl` | job body | Builds the job from `LaBrRun`. Never edited |
| `fcl/LaBrPipeline.fcl` | run file | The file you copy and edit |
| `fcl/LaBrPipeline_PoissonSpectrumPlusIntrinsic.fcl`, `fcl/LaBrPipeline_UpstreamDataset.fcl` | run files | Worked examples |
| `fcl/DumpSTMParticles.fcl` | job | Runs `STMParticleDump` on upstream art files |
| `src/makeLinesFcl.C` | ROOT macro | Turns a histogram with one bin per line (e.g. `stmdet/hLinesCut`) into a run file |
| `src/compareLaBrSpectra.C` | ROOT macro | Old-vs-new (or any A-vs-B) comparison |
| `src/plotLaBrBinary.py` | Python | Plots / exports the raw `.bin` stream, finds pulses, overlays truth |

### Existing pieces the pipeline reuses (unchanged)
| File | Role in the pipeline |
|---|---|
| `src/LaBrIntrinsicSource_module.cc` + `inc/LaBrIntrinsicDecays.hh` | La-138 and Ac-227 chain decay products, emitted inside the crystal |
| `src/GenParticleMerger_module.cc` | Concatenates generator outputs for g4run |
| `src/LaBrWaveformsFromStepPointMCsTest_module.cc` | The LaBr digitizer |
| `src/LaBrTruthPulses_module.cc`, `src/TagPileupPulses_module.cc` | Truth pulse per decay; pileup tagging |
| `src/LaBrTruthSpectrum_module.cc` | Truth energy per event and per decay |
| `inc/LaBrTruthHelpers.hh` | Shared alpha quenching, primary tracing, source classification, decay grouping |
| `src/BinaryWaveformWriter_module.cc` | Raw int16 ADC stream |
| `fcl/prolog.fcl` | `LaBrDigitization`, `LaBrPileup`, `STMDAQParameters`, `ComponentPositions`, `ShiftVD101Steps`, `LaBrIntrinsicActivity` |
| `Offline/Mu2eG4/fcl/prolog.fcl` | The default `g4run` configuration |
| `Offline/Mu2eG4/geom/STM_v09.txt` (via `geom_common.txt` -> `geom_run1_a.txt`) | The Geant4 STM geometry |

---

## 3. How a job is assembled (the FHiCL layer)

### 3.1 The three-part run file
Every run file has the same shape:

```fcl
#include "Offline/STMMC/fcl/LaBrPipelineProlog.fcl"   # (1) defines LaBrPipeline.* and default LaBrRun
BEGIN_PROLOG
LaBrRun.nEvents : 10000                             # (2) your overrides, applied in order
LaBrRun.primary : @local::LaBrPipeline.Primary.gun
...
END_PROLOG
#include "Offline/STMMC/fcl/LaBrPipelineJob.fcl"      # (3) builds the job from the final LaBrRun
```

1. `LaBrPipelineProlog.fcl` includes the standard services and producers (which define
   `g4run`) and the STMMC `prolog.fcl`. It then defines the building blocks one piece at a time
   (`LaBrPipeline.Beam`, `.Emission`, `.Primary`, `.intrinsic`, `.DigiLaBr`, `.TruthPulses`,
   `.TagPileup`, `.TruthLaBr`, `.BinaryWriter`). Each piece is defined as its own statement so
   later pieces can `@local`-reference earlier ones: FHiCL cannot reference a key inside a
   table that is still being defined. Finally it defines the default `LaBrRun`.
2. Your block changes `LaBrRun`. FHiCL applies these statements **in order**. So replace a whole
   table (`LaBrRun.primary : ...`) **before** changing keys inside it
   (`LaBrRun.primary.E : ...`); the other order silently throws your key away.
3. `LaBrPipelineJob.fcl` contains no includes of its own. It reads the final `LaBrRun` and
   `LaBrPipeline.*` values and builds `process_name`, `source`, `services`, `physics` and
   `outputs`.

`LaBrPipelineJob.fcl` must be included last, and must not be included together with a second
copy of the prolog. Re-including the prolog would reset `LaBrRun` to its defaults.

### 3.2 How `LaBrRun` maps onto the job
| `LaBrRun` key | Used in `LaBrPipelineJob.fcl` as |
|---|---|
| `processName` | `process_name` (no underscores: art forbids them) |
| `nEvents` | `source.maxEvents` **and** `physics.producers.TagPileup.nMerge` |
| `baseSeed` | `services.SeedService.baseSeed` |
| `primary`, `background`, `intrinsic` | the producer tables with those labels |
| `sources` | `genCombined.inputTags` **and** the start of `simulation_path` |
| `truthProducers` | the end of `simulation_path` (default `[TruthPulses, TagPileup]`; `[]` skips them) |
| `endModules` | `end_path` (any of `TruthLaBr`, `BinaryWriter`, `fullOutput`) |
| `artFile`, `histFile`, `binaryFile` | `outputs.fullOutput.fileName`, `services.TFileService.fileName`, `BinaryWriter.fileName` |

The path is
`simulation_path : [ <sources...>, genCombined, g4run, DigiLaBr, <truthProducers...> ]`.
The three generator slots are always defined as producers, but art only builds modules that are
on a path. So a slot not listed in `sources` is never constructed, and its settings (e.g. a
dataset file that does not exist) do not matter.

### 3.3 Two ways to change a building block
- **In your prolog block, before the job include** (preferred):
  `LaBrPipeline.DigiLaBr.excessNoiseFactor : 1.5`. This works because the job copies
  `@local::LaBrPipeline.DigiLaBr` after your block has run.
- **After the job include**: `physics.producers.DigiLaBr.excessNoiseFactor : 1.5` at the very end
  of the run file. This works for anything, including `g4run` and `services`.

**Values are copied when referenced, not linked.** `LaBrPipeline.DigiLaBr` took its values from
`LaBrDigitization.*` when the prolog file was read. Setting `LaBrDigitization.excessNoiseFactor`
in your block afterwards does **nothing** to the job. Always override the `LaBrPipeline.*` (or
`physics.*`) key itself.

---

## 4. One micropulse, step by step

### 4.1 The event
- The source is `EmptyEvent`. Each art event is one micropulse window of **1695 ns**, and
  art event numbers start at **1**.
- Nothing carries time from one event to the next except two things:
  - the digitizer's charge buffer (4.5.6);
  - the truth timeline (4.6).
- Services: `Services.Core` + `Services.SimOnly`. The GeometryService uses
  `Offline/Mu2eG4/geom/geom_common.txt`. SeedService is `autoIncrement` from `LaBrRun.baseSeed`:
  each module asking for a seed gets the next number, in the order the modules are built.

### 4.2 Generator: `STMParticleSource` (labels `primary`, `background`)

**At construction (once per job):**
1. **`mode`**:
   - `"gun"` needs `E` (kinetic energy, MeV).
   - `"spectrum"` builds a normalised cumulative distribution (CDF):
     - `lines`: CDF over the lines (discrete);
     - `histogram`: CDF over bins given as N+1 `binEdges` and N `binContents` (counts per bin,
       not per MeV);
     - `rootHist`: the same built from a TH1's bins, x axis in MeV.
     All weights must be >= 0 with a positive total.
   - `"dataset"` loads the whole TTree into memory (4.2.1).
2. **`beam`** (gun/spectrum only):
   - start point `(x, y, z)` mm;
   - optional uniform disk `radius` in the x-y plane;
   - unit `direction` (normalised);
   - optional cone `coneHalfAngle` in degrees (0 = pencil beam, 180 = isotropic).
3. **`emission`**:
   - `poisson : true`: mean = `nPerEvent`, or `rateHz` x `windowTime` x 1e-9 if `rateHz` is given.
   - `poisson : false`: exactly round(`nPerEvent`) particles per event; `rateHz` is then an error.
   - `timing`:
     - `"uniform"`: U[0, `windowTime`).
     - `"fixed"`: `time`.
     - `"dataset"`: the upstream time folded into the window, `fmod(t, windowTime)`; dataset
       mode only.
4. The random engine is created from SeedService (`CLHEP::RandFlat` + `RandPoissonQ`).

**Each event:**
1. **How many.** N = Poisson(mean), or the fixed count.
2. **Gun / spectrum, for each of the N particles:**
   - **Energy.** `gun`: E. `lines`: pick a line from the CDF. `histogram`/`rootHist`: pick a bin
     from the CDF, then a uniform energy within that bin.
   - **Position.** The start point, plus a point uniform in the disk when `radius` > 0
     (r = R sqrt(u), phi uniform).
   - **Direction.** The axis, or uniform in solid angle inside the cone:
     cos(theta) uniform in [cos(theta_max), 1], phi uniform.
   - **Momentum.** |p| = sqrt(Ek (Ek + 2m)), with m from the particle data table
     (GlobalConstantsService). E_total = sqrt(p^2 + m^2).
   - **Time.** From `timing`.
   - **Output.** `GenParticle(pdgId, GenId::particleGun, pos, p4, t)`.
3. **Dataset, for each of the N sampling units:**
   - pick a unit uniformly at random, with replacement (a particle, or a whole upstream event);
   - pick a random start time s in [0, window);
   - emit each particle of the unit with its stored pdgId, position and momentum;
   - time: with `timing : "uniform"` it is t = s + (t_i - t_earliest in the unit), so particles
     keep their relative timing. Particles later than 2 x 1695 ns are dropped by the digitizer
     (4.5.1). With `timing : "dataset"` it is `fmod(t_i, window)`.

`GenId::particleGun` makes every particle count as source 0 ("external") in the truth modules.

#### 4.2.1 Dataset loading and the transform
- **Branches read.** `pdgId` (int), `x y z px py pz time` (double), and for `sampleUnit : "event"`
  also `eventId` (Long64). The default tree is `DumpParticles/particles`, as written by
  `STMParticleDump`.
- **Filters**, applied at load:
  - `pdgIds` (keep list);
  - kinetic energy in [`EMin`, `EMax`] (`EMax` <= 0 means no upper limit);
  - entries whose PDG code isn't in the particle table are skipped and counted in the log.
- **Grouping.** With `sampleUnit : "event"`, consecutive entries with the same `eventId` form one
  unit (the tree must be ordered by event, as `STMParticleDump` writes it).
- **Optional `transform`** (maps an upstream plane onto the aperture), for each particle:
  - d = (pos - inputCentre) with its z set to 0. The particle is dropped if |d| > `inputRadius`
    (when > 0).
  - pos' = outputCentre + s d (so z' = outputCentre.z), with s = `radialScale`.
  - p'_x = s p_x, p'_y = s p_y, and p'_z = sign(p_z) sqrt(|p|^2 - p'_x^2 - p'_y^2): |p| is kept.
  - This is the same idea as `ShiftVirtualDetectorStepPointMCs` in the production resampler, but
    without its special handling of the HPGe absorber shadow. The default maps VD101 (r = 200 mm)
    onto the LaBr SSC aperture (r = 3.989 mm): s = 3.989/200.

#### 4.2.2 `STMParticleDump` (making a dataset)
- For each input event, it keeps the StepPointMCs in the given collection whose
  `virtualDetectorId` is in `virtualDetectorIds` (default `[101]`), filtered by `pdgIds`.
- With `forwardOnly`, only p_z > 0 is kept.
- With `firstCrossingOnly`, only the earliest crossing per (track, VD) is kept, which drops
  re-crossings after backscatter.
- Tree `particles`: eventId (sequential), run/subRun/event, virtualdetectorId, pdgId, trackId,
  x y z [mm], px py pz [MeV/c], time [ns], E [MeV kinetic].
- Tree `info`: nEvents, the number of input events read, for normalisation.

### 4.3 Generator: `LaBrIntrinsicSource` (label `intrinsic`)
- **Crystal volume**: V = pi r^2 (2 h) from `radius`, `halfLength` (default 347.5 cm^3).
- **La-138 specific activity**: computed from the composition (density 5.08 g/cm^3, molar mass
  378.62, 0.95 La per formula unit, abundance 8.881e-4, T1/2 1.03e11 yr), giving about
  1.454 Bq/cm^3 (505 Bq). `la138SpecificActivity` overrides it.
- **Ac-227 chain**: 8 members, each at `ac227ChainActivity` x its branching fraction (secular
  equilibrium), or `memberSpecificActivities` (8 values).
- **Mean decays per event** for each species: A_specific x V x `activityScale` x `windowTime` x
  1e-9. At `activityScale` 1 that is 8.6e-4 La-138 decays per micropulse; at 1000, 0.86.
- **Each decay** is at one vertex uniform in the cylinder and one time uniform in the window. It
  emits its decay products directly (Geant4 radioactive decay is off):
  - La-138 EC (65.5%): 1435.8 keV gamma + Ba atomic relaxation (K/L/M vacancy: X-rays, or an
    Auger electron for the rest);
  - La-138 beta- : 788.7 keV gamma + a beta electron (263.7 keV endpoint,
    unique-2nd-forbidden shape);
  - Ac-227 chain members: alphas / betas / gammas (see `inc/LaBrIntrinsicDecays.hh`).
- **Tags**: products are tagged `GenId::STMLa138` or `GenId::STMAc227`.
- `includeLa138` / `includeAc227` switch each part off.

### 4.4 Merge and Geant4
- **`genCombined`** (`GenParticleMerger`) concatenates the `GenParticleCollection`s of the labels
  in `LaBrRun.sources`, in list order.
- **`g4run`** (Mu2eG4) is `@local::g4run` from `Offline/Mu2eG4/fcl/prolog.fcl`, with these
  overrides from the job file:
  - `physics.physicsListName : "QGSP_BERT_EMZ"` (EM option 4, accurate for keV-MeV photons);
  - `physics.turnOnRadioactiveDecay : false` (the generators already emit the decay products);
  - `SDConfig.enableSD : [ STMDet ]` (only the STM crystal sensitive detector);
  - `inputs.primaryTag : "genCombined"` (`primaryType` is the default, GenParticles).
  - The defaults kept include `minRangeCut` 0.010 mm and no Mu2eG4 common cut.
- **Geometry**: the full Mu2e world from GeometryService (Section 7). Particles are transported
  through everything: SSC collimator, shielding, the LaBr endcap/window, the crystal.
- **Output**: `g4run:STMDet`, one StepPointMC per step in **either** STM crystal (HPGe and LaBr
  share the collection), plus the SimParticle tree (parents, start vertex/time, GenParticle
  link) that the truth modules use.
- **Times**: `step.time()` is the global time in ns, relative to the event start: the
  GenParticle time plus the flight time.

### 4.5 Digitizer: `LaBrWaveformsFromStepPointMCsTest` (label `DigiLaBr`)
Units inside the digitizer are **photoelectrons (pe)** until the final conversion to ADC.

**Constants built at construction** (with the pipeline values):
- tADC = 1000/fADC = 1000/370.37 = **2.700 ns**.
- nADCs = int(1695/tADC) + 1 = **628** samples per micropulse.
- Buffer `_chargeCollected` = nADCs x `microspillBufferLengthCount` (2) = 1256 bins.
- pe/MeV = `lightYield` x `quantumEfficiency` = 63000 x 0.35 = **22050**.
- Pulse CDF table: the single-photoelectron shape s(t) = (1 - e^(-t/tau1))^p e^(-t/tau2)
  (tau1 = 22.3543, tau2 = 28.0007 ns, p = 6.6596), integrated with the trapezoidal rule on a
  0.05 ns grid out to 40 tau2 and normalised to 1. It is tabulated because it has no closed form
  for non-integer p.
- **Peak fraction**: the largest share of one deposit that falls in a single bin,
  max_k [CDF((k+1) tADC) - CDF(k tADC)] = **0.04114** at 370.37 MHz.
- **Gain** (peak calibration, `peakKeVPerADC` = 0.57):
  chargeToADC = (1000 / pe_per_MeV) / (peakKeVPerADC x peakFraction) = **1.934 ADC per pe**.
  So a deposit of E keV has a peak of E/0.57 ADC: **1809 keV -> 3174 ADC**, 347 keV -> 609 ADC.
  (With `peakKeVPerADC` unset it uses the integral calibration `EnergyPerADCBin` instead.)

**Each event:**
1. **Step selection**, for every StepPointMC with ionizing energy > 0:
   - reject x < -3925 mm (the HPGe side);
   - reject step.time() > 2 x 1695 ns (beyond the buffer);
   - change to crystal-local coordinates: subtract the hard-coded centre
     (-3863.36, 0, 40891.50), then shift so the front face is at z = 0;
   - reject if outside the cylinder (R <= 38.1, 0 <= z <= 76.2, with 0.1 mm tolerance).
   The first 20 rejects are printed with their world coordinates, and every reject is counted in
   the end-of-job "loss breakdown".
2. **Light**: E_light = ionizing E, times the alpha quenching factor for alpha steps
   (`alphaQuenchingFactor` + slope x (E0 - 6 MeV), E0 = the alpha's initial energy).
3. **Photostatistics**: N_pe = Poisson(E_light x 22050 / MeV).
4. **PMT excess noise**: signal = N_pe + Gauss(0, 1) x sqrt((F - 1) N_pe), F = 1.27. Together
   this gives sigma/E = sqrt(F / N_pe):
   - 662 keV: 2.20% FWHM (14.5 keV);
   - 1809 keV: sigma 10.2 keV (24.0 keV FWHM).
5. **Pulse placement**: start bin t0 = floor((step.time() + timeOffset) / tADC). Bin t0 + k
   receives -signal x [CDF((k+1) tADC) - CDF(k tADC)]; the pulse is negative-going. This stops at
   the buffer end or when less than 0.5 pe is left. The peak is about 50 ns (about 18 samples)
   after the start.
6. **Carry-over**: every step deposits into the 2-micropulse buffer. After digitizing, bins
   628-1255 become bins 0-627 of the next event. So a deposit late in micropulse N continues
   smoothly into N+1, and a particle arriving at 1695-3390 ns lands in the next micropulse.
7. **`decayCharge`**: a straight copy (direct PMT readout, no preamplifier).
8. **Noise**: Gaussian with `NoiseSD`, off (0) in the pipeline. See the warning in
   [Section 10](#10-validation-status-and-known-limitations) before turning it on.
9. **ADC**: ADC_i = round(charge_i x chargeToADC), floored at -32767 (saturation at about
   18.7 MeV equivalent peak). Stored as int16.
10. **Output**:
    - one `STMWaveformDigi` (628 samples) as `DigiLaBr:LaBr`;
    - with `makeTTree`, 628 rows in `DigiLaBr/ttree`: `eventId`, `time` (sample counter),
      `ADC`, `chargeCollected`, `chargeDecayed` (pe).

The digitizer's random engine is a `std::default_random_engine` that is **never seeded**: it
replays the same sequence in every job, whatever `baseSeed` is. This is a known issue, noted and
deliberately not yet fixed.

### 4.6 Truth: `LaBrTruthPulses` -> `TagPileup`, and `TruthLaBr`
All three use the **same step selection** as the digitizer (x cut + the same cylinder), and they
share the helpers in `inc/LaBrTruthHelpers.hh`:
- **Primary**: walk `SimParticle::parent()` to the particle with no parent. `originParticle()`
  is **not** used, because it returns ordinary secondaries themselves.
- **Source**: `GenId::STMAc227` -> 2, `GenId::STMLa138` -> 1, otherwise vertex inside the
  crystal -> 1, otherwise -> 0 (external).
- **Decay key**: (source, primary start vertex rounded to 1 um, start time rounded to 1 ps). All
  products of one decay share it, so they form **one** pulse. A La-138 EC gamma + X-ray + Auger
  electron is one 1468 keV pulse, not three. Two gun photons with the same vertex and time
  (possible with `single`, `nPerEvent` > 1, a point beam and fixed timing) also merge, which
  matches what the detector would see.
- **Energies** are light-equivalent (the same alpha quenching as the digitizer).

**`LaBrTruthPulses`**, per event:
- one `STMMCPulse` per decay key: energy = summed deposit (keV), time = the earliest deposit;
- the time goes on an **absolute** axis: t_abs = (events so far) x 1695 ns + local time;
- `isPileup` is left false.

**`TagPileup`** (`TagPileupPulses`):
- collects the pulses of `nMerge` events (= `LaBrRun.nEvents`, so the whole run is one timeline)
  and sorts them by time;
- for each pulse, w_after(E) = tau2 x ln(peakADCPerKeV x E / thresholdADC), floored at
  `minWindow` (20 ns):

  | Energy | w_after |
  |---|---|
  | 66 keV | 88 ns |
  | 347 keV | 134 ns |
  | 844 keV | 159 ns |
  | 1809 keV | 181 ns |

- two pulses dt apart are both tagged pileup if dt < w_after(earlier) or dt < `wBefore` (20 ns);
- writes the tagged collection, and the `TagPileup/pulses` tree (time, energy, source, eventId,
  primaryId, isPileup, nOverlap, dtNearest, window);
- if the job ends mid-batch, the leftover pulses are still tagged and written to the tree at
  end of job.

**`TruthLaBr`** (`LaBrTruthSpectrum`, analyzer), writes:
- `TruthLaBr/ttree`: one entry per event (edepTotal, edepGun, edepLa138, edepAc227, keV);
- `TruthLaBr/decays`: one entry per decay (energy, source, local time of the first deposit).

### 4.7 Writing
- **`BinaryWriter`** appends each event's 628 int16 samples to `LaBrRun.binaryFile`: headerless,
  little-endian, back to back. Size = 1256 bytes x events.
- **`fullOutput`** writes every product (generators, g4run, digis, truth) to `LaBrRun.artFile`.
- **TFileService** writes the trees above to `LaBrRun.histFile`.

---

## 5. The numbers that tie the stages together

| Quantity | Value | Where it comes from |
|---|---|---|
| Micropulse | 1695 ns | hard-coded in the digitizer; `LaBrPileup.micropulseTime`; `emission.windowTime` |
| Sample spacing | 2.700 ns | 1000 / `STMDAQParameters.samplingFrequencies.LaBr` (370.37 MHz) |
| Samples per micropulse | 628 | int(1695/2.7) + 1 |
| Waveform time per event | 1695.6 ns | 628 x 2.7. **Not** 1695: map truth to samples per event, never with a single global time -> sample formula (it drifts by about 2 samples per 10k events) |
| Photoelectrons | 22050 per MeV | 63000 photons/MeV x QE 0.35 |
| Peak fraction | 0.04114 | computed from the pulse shape and fADC |
| Gain | 1.934 ADC per pe; 0.57 keV per ADC of peak | `peakKeVPerADC` |
| 1809 / 347 / 66 keV peak | 3174 / 609 / 116 ADC | E / 0.57 |
| Saturation | -32767 ADC (about 18.7 MeV) | int16 |
| Resolution | 2.20% FWHM at 662; sigma 10.2 keV at 1809 | sqrt(F/N_pe), F = 1.27 |
| Pulse peak delay | about 50 ns after the start | pulse shape |
| 20 kHz | 0.0339 per micropulse | rateHz x 1695e-9 |
| La-138 | 505 Bq, 8.6e-4 per micropulse (x `activityScale`) | composition and volume |
| Truth event origin | (eventId - 1) x 1695 ns | `LaBrTruthPulses` |

---

## 6. How to use it: recipes

Start from a copy of `fcl/LaBrPipeline.fcl`. Put it in your run directory and run it with
`mu2e -c $PWD/MyRun.fcl`. Every snippet below goes in the `BEGIN_PROLOG` block unless it says
"after the include". Keep fcl files ASCII-only.

### 6.1 Run length, names, seed
```fcl
LaBrRun.processName : "MyRun"                 # letters/digits only
LaBrRun.nEvents     : 100000                  # micropulses (100000 = 169.5 ms)
LaBrRun.baseSeed    : 9                       # changes generators + Geant4 (not digitizer noise)
LaBrRun.artFile     : "sim.MyRun.art"
LaBrRun.histFile    : "nts.MyRun.root"
LaBrRun.binaryFile  : "labr_MyRun_ADC.bin"
LaBrRun.endModules  : [ TruthLaBr, BinaryWriter, fullOutput ]   # drop fullOutput to skip the big .art
```
Change the event count with `LaBrRun.nEvents`, not `mu2e -n`, so the pileup batch still covers
the whole run.

### 6.2 Photon gun
```fcl
LaBrRun.primary : @local::LaBrPipeline.Primary.gun
LaBrRun.primary.E : 0.662                     # MeV kinetic
LaBrRun.primary.pdgId : 22                    # any PDG code in the particle table, e.g. 2112 neutron, 11 electron
LaBrRun.primary.emission : @local::LaBrPipeline.Emission.single   # 1 per event at t = 0
```

### 6.3 Spectrum
```fcl
LaBrRun.primary : @local::LaBrPipeline.Primary.spectrum
# lines
LaBrRun.primary.spectrum : { type : "lines"  lineEnergies : [0.347, 1.809]  lineIntensities : [3, 1] }
# or a histogram typed in (N+1 edges, N contents = counts per bin)
LaBrRun.primary.spectrum : { type : "histogram"  binEdges : [0.0, 0.5, 1.0, 2.0]  binContents : [10, 5, 1] }
# or a TH1 from a file (x axis in MeV; flat within each bin)
LaBrRun.primary.spectrum : { type : "rootHist"  fileName : "/abs/path/spec.root"  histName : "dir/hE" }
```
- `fileName` can be a path (absolute, or relative to where you run) or a name on
  `MU2E_SEARCH_PATH`.
- For a histogram whose bins **are** lines (one bin per line, labelled in keV, like
  `stmdet/hLinesCut`), use `makeLinesFcl.C`, not `rootHist` (see 6.9). `rootHist` would spread
  each line uniformly across its bin.

### 6.4 Poisson and rate
```fcl
LaBrRun.primary.emission : @local::LaBrPipeline.Emission.poisson
LaBrRun.primary.emission.rateHz : 20000.0     # mean per event = rate x 1695 ns
# or give the mean directly (and remove rateHz, which takes priority):
#   LaBrRun.primary.emission : { poisson : true  nPerEvent : 0.05  timing : "uniform"  windowTime : 1695.0 }
```
Keep `windowTime` at 1695 (the digitizer's micropulse is fixed).

### 6.5 Beam position, size and direction
```fcl
LaBrRun.primary.beam.x : -3863.4              # mm, Mu2e coordinates
LaBrRun.primary.beam.z : 40404
LaBrRun.primary.beam.radius : 3.5             # uniform disk; the LaBr SSC aperture is r = 3.989 mm
LaBrRun.primary.beam.direction : [0, 0, 1]
LaBrRun.primary.beam.coneHalfAngle : 2.0      # deg; 180 = isotropic point source
```
The disk lies in the x-y plane whatever the direction. The default start (z = 40404, on the LaBr
aperture axis) is upstream of the SSC, so Geant4 transports the photons through the real
collimator.

### 6.6 Self-scintillation
```fcl
LaBrRun.sources : [ primary, intrinsic ]
LaBrRun.intrinsic.includeAc227 : false        # La-138 only
LaBrRun.intrinsic.activityScale : 100.0       # boosts statistics AND pileup: rates are unphysical unless 1
LaBrRun.intrinsic.la138SpecificActivity : 1.425   # literature value instead of the computed one (Bq/cm^3)
LaBrRun.intrinsic.ac227ChainActivity : 0.03       # Bq/cm^3 per chain member
```
Only the intrinsic source: `LaBrRun.sources : [ intrinsic ]`. The alpha quenching (0.30) is
already applied in the digitizer and in both truth modules.

### 6.7 Upstream dataset
1. Make the tree from upstream art files, e.g. the BeamToVD101 / NeutralsToVD101 outputs:
   ```bash
   mu2e -c Offline/STMMC/fcl/DumpSTMParticles.fcl -S upstream_files.txt   # -> nts.STMParticleDump.root
   ```
   In that fcl you can change `StepPointMCsTag`, `virtualDetectorIds`, `pdgIds`, `forwardOnly`
   and `firstCrossingOnly`.
2. Use it:
   ```fcl
   LaBrRun.primary : @local::LaBrPipeline.Primary.dataset
   LaBrRun.primary.dataset.fileName : "/abs/path/nts.STMParticleDump.root"
   LaBrRun.primary.dataset.sampleUnit : "event"      # keep each upstream event's particles together
   LaBrRun.primary.dataset.pdgIds : [ 2112 ]         # optional: neutrons only
   LaBrRun.primary.dataset.EMin : 0.01               # MeV
   LaBrRun.primary.emission.nPerEvent : 0.8          # SET from the normalisation below
   ```
3. **Normalisation** is not automatic:
   - `sampleUnit "event"`: nPerEvent = upstream events per micropulse, e.g.
     (POT per micropulse) / (POT per upstream event).
   - `sampleUnit "particle"`: that, times `particles->GetEntries() / info.nEvents`.
4. **Transform**: the default maps VD101 onto the SSC aperture. If your tree was recorded at
   another plane, set `transform.inputCentre`, `inputRadius` and `radialScale`, or remove the
   transform to use the positions as recorded:
   ```fcl
   LaBrRun.primary.dataset : { fileName : "..."  treeName : "DumpParticles/particles"  sampleUnit : "event"  EMin : 0  EMax : 0 }
   ```
5. A tree from another program works if it has the branches in 4.2.1.

### 6.8 Signal + background + intrinsic in one job
```fcl
LaBrRun.primary : @local::LaBrPipeline.Primary.gun
LaBrRun.primary.emission : @local::LaBrPipeline.Emission.poisson
LaBrRun.background.dataset.fileName : "/abs/path/nts.STMParticleDump.root"
LaBrRun.background.emission.nPerEvent : 0.8
LaBrRun.sources : [ primary, background, intrinsic ]
```
Both STMParticleSource slots tag their particles as source 0, so truth cannot tell `primary` from
`background` apart (see Section 10).

A third STMParticleSource needs its own label, added after the include:
```fcl
# after #include "Offline/STMMC/fcl/LaBrPipelineJob.fcl":
physics.producers.line844 : @local::LaBrPipeline.Primary.gun
physics.producers.line844.E : 0.844
physics.producers.line844.emission : @local::LaBrPipeline.Emission.poisson
physics.producers.line844.emission.rateHz : 500
physics.producers.genCombined.inputTags : [ primary, line844 ]
physics.simulation_path : [ primary, line844, genCombined, g4run, DigiLaBr, TruthPulses, TagPileup ]
```
`@local::` can still be used here, because the prolog names remain available until the end of
the file.

### 6.9 Line mix from a histogram (e.g. the SSC photon study)
```bash
root -l -b -q 'Offline/STMMC/src/makeLinesFcl.C("gun_check_FOILOFF_w.root","stmdet/hLinesCut","FOILOFF_20kHz",20000)'
mu2e -c $PWD/LaBrSSCLines_FOILOFF_20kHz.fcl
```
- Arguments: input file, histogram (one bin per line, bin labels = energy in keV), tag, rateHz
  (0 = one photon per event at t = 0).
- The generated fcl writes the `.bin` too. Edit it like any run file.

### 6.10 Detector response
Override the `LaBrPipeline.DigiLaBr` keys (see 3.3):
```fcl
LaBrPipeline.DigiLaBr.lightYield : 73000          # photons/MeV
LaBrPipeline.DigiLaBr.quantumEfficiency : 0.35
LaBrPipeline.DigiLaBr.excessNoiseFactor : 1.27    # resolution: sigma/E = sqrt(F / (E x lightYield x QE))
LaBrPipeline.DigiLaBr.pulseTau1 : 22.3543         # ns   pulse shape (1-e^-t/tau1)^p e^-t/tau2
LaBrPipeline.DigiLaBr.pulseTau2 : 28.0007         # ns
LaBrPipeline.DigiLaBr.pulseRiseSharpening : 6.6596
LaBrPipeline.DigiLaBr.peakKeVPerADC : 0.57        # gain: keV per ADC of peak height
LaBrPipeline.DigiLaBr.fADC : 300.0                # MHz: changes samples per event (1695/3.333 + 1 = 509)
LaBrPipeline.DigiLaBr.microspillBufferLengthCount : 3   # accept arrivals up to 3 x 1695 ns
LaBrPipeline.DigiLaBr.alphaQuenchingFactor : 0.30
```
Some of these have copies elsewhere that must change with them (Section 9): the pulse decay,
gain, fADC, micropulse and alpha quenching.

### 6.11 Pileup tagging
```fcl
LaBrPipeline.TagPileup.thresholdADC : 10.0        # set to your pulse-finder threshold
LaBrPipeline.TagPileup.wBefore : 20.0             # ns, rise-side window
LaBrPipeline.TagPileup.minWindow : 20.0           # ns
LaBrRun.truthProducers : []                       # skip truth pulses + pileup entirely
```
To re-tag an existing art file without Geant4, use `LaBrTruthPulsesFromArt.fcl`. Its
`process_name` must differ from the input's, and `nMerge` must equal the input event count.

### 6.12 Geant4 settings (after the include)
```fcl
physics.producers.g4run.physics.physicsListName : "QGSP_BERT_EMZ"
physics.producers.g4run.physics.minRangeCut : 0.010      # mm
physics.producers.g4run.SDConfig.enableSD : [ STMDet, virtualdetector ]   # also record VD steps
```

---

## 7. Changing the geometry

Geometry lives in **two** places, and a change must be made in both:
- the Geant4 world, built by GeometryService from text files; and
- **hard-coded copies** of the LaBr crystal position and size in C++ and fcl, which the digitizer
  and truth use to select steps and which the intrinsic source uses to place decays.

### 7.1 Where the Geant4 geometry comes from
`services.GeometryService.inputFile` (`Services.Core`) is `Offline/Mu2eG4/geom/geom_common.txt`.
That includes `geom_run1_a.txt`, which includes `STM_v09.txt`. The LaBr-relevant parameters
there:

| Parameter (`STM_v09.txt`) | Value | Meaning |
|---|---|---|
| `stm.LaBr.build` | true | build the LaBr detector |
| `stm.LaBr.crystalMaterial` / `windowMaterial` / `wallMaterial` | LaBr3Ce / G4_Al / G4_Al | materials |
| `stm.LaBr.CrystalR`, `stm.LaBr.CrystalL` | 38.1, 76.2 mm | crystal radius and length (3" x 3") |
| `stm.LaBr.EndcapR`, `EndcapL`, `WindowD`, `EndcapD`, `AirD` | 48.0, 66.3, 1.5, 1.5, 4.0 mm | housing |
| `stm.LaBr.Z_LaBr` | 330 mm | detector position downstream of the SSC reference |
| `stm.LaBr.offset_LaBr` | 0 mm | extra x offset of the LaBr |
| `stm.STM_SSC.offset_Spot` | 40.64 mm | x offset of each detector's spot from the beam axis |
| `stm.STM_SSC.FrontToWall` | 1049.2 mm | SSC front face, upstream of the east wall |
| `stm.STM_SSC.Aperture_LaBr1`, `Aperture_LaBr2` | 50, 150 | LaBr SSC apertures. 50 appears to be an **area** in mm^2: sqrt(50/pi) = 3.989 mm, exactly the radius the production shift uses |
| `stm.STM_SSC.material`, depths | G4_W, 50.8 mm | collimator |
| `stm.STM_Absorber.build` | true | HPGe absorber |

In `GeometryService/src/STMMaker.cc` the LaBr is placed at
`STMShieldingRef + (offset_Spot + offset_LaBr, 0, Z_LaBr)`, with
`STMShieldingRef = (beam axis at the east wall) - (0, 0, FrontToWall)`. So the LaBr sits at
x = -3904 + 40.64 = **-3863.36**, and moving `Z_LaBr` moves it along z.

### 7.2 Overriding Geant4 geometry parameters
GeometryService allows a later definition to replace an earlier one (`allowReplacement` is true).
So write a small geometry file that includes the standard one and redefines what you want:

```
// Offline/STMMC/geom/labr_study.txt
#include "Offline/Mu2eG4/geom/geom_common.txt"
double stm.LaBr.Z_LaBr    = 400;      // move the LaBr 70 mm further downstream
double stm.LaBr.CrystalR  = 25.4;     // 2" x 2" crystal
double stm.LaBr.CrystalL  = 50.8;
```

Then point the job at it, after the include in your run file:
```fcl
services.GeometryService.inputFile : "Offline/STMMC/geom/labr_study.txt"
```
Put the file somewhere under your Muse work area, as with all `Offline/...` paths, so the
GeometryService file lookup finds it. Geometry text files are read at run time: no rebuild is
needed for this part.

### 7.3 The hard-coded copies you must update
A Geant4 change alone is **not enough**. These copies do not read GeometryService:

| What | File : line | Used for |
|---|---|---|
| Crystal centre x, y, z | `src/LaBrWaveformsFromStepPointMCsTest_module.cc:140-142` | digitizer step selection |
| Crystal L, R | same file `:145-146` | digitizer bounds |
| HPGe/LaBr x cut (-3925) | same file `:400` | digitizer rejects HPGe steps |
| Centre, L, R, x cut | `src/LaBrTruthPulses_module.cc:81-86` | truth pulses |
| Centre, L, R, x cut | `src/LaBrTruthSpectrum_module.cc:77-85` | truth spectrum |
| `LaBrPipeline.intrinsic` x, y, z, radius, halfLength | `fcl/LaBrPipelineProlog.fcl:122-126` | where La-138 / Ac-227 decays are placed (and the volume, so the activity) |
| `LaBrPipeline.Beam.aperture` x, z | `fcl/LaBrPipelineProlog.fcl:33-35` | gun start point |
| Dataset `transform` (outputCentre, inputRadius, radialScale) | `fcl/LaBrPipelineProlog.fcl:105-107` | dataset mapping onto the aperture |
| `ComponentPositions.LaBr`, `ShiftVD101Steps.LaBrUpStr`, `OutputRadius` | `fcl/prolog.fcl:34-35, 57-65` | aperture position/radius used by the above |

The C++ changes need `muse build -j8`, and you must check that the build **linked** (a failed
link leaves the old plugins in place). The fcl changes can go in your run file instead of
editing the prolog:
```fcl
LaBrRun.intrinsic.z : 40960.0
LaBrRun.intrinsic.radius : 25.4
LaBrRun.intrinsic.halfLength : 25.4
LaBrRun.primary.beam.z : 40404
```

**The 101.6 mm offset.** The digitizer's crystal z (40891.50) is **not** what `STMGeometryDump`
reports (40993.1). It was corrected empirically by -101.6 mm to match where Geant4 actually puts
the steps, and the cause is unknown. When you move the detector, do not copy the dump value:
find the new centre as in 7.4.

### 7.4 Checking a geometry change
1. `mu2e -c Offline/STMMC/fcl/STMGeometryDump.fcl` prints the positions GeometryService builds.
   Expect the same 101.6 mm z discrepancy.
2. Run a short gun job (`LaBrRun.nEvents : 1000`, 1809 keV, `single`). At end of job the
   digitizer prints its **loss breakdown**. `rejected: out-of-bounds` should be about 0. If not,
   the first 20 rejects are printed with world coordinates. Their spread shows the real crystal
   centre (the midpoint of the step z range) and radius.
3. With `LaBrPipeline.TruthLaBr.makeStepTree : true`, `TruthLaBr/steps` holds the x, y, z of
   every kept deposit, which you can histogram to check the crystal volume directly.
4. Update the constants in 7.3, rebuild, and repeat until nothing is rejected and the photopeak
   efficiency is sensible.

A smaller crystal also changes the La-138 activity, because `LaBrIntrinsicSource` computes it
from its own `radius` and `halfLength`. So set those to the new size too.

---

## 8. Outputs and the tools that read them

| Output | Contents | Read with |
|---|---|---|
| `LaBrRun.binaryFile` (.bin) | int16 little-endian, 628 samples/event, back to back | `src/plotLaBrBinary.py` (overview, windows, `--plain`, `--pulses`, `--export`, `--spectrum`, `--truth`); `np.fromfile(f, dtype="<i2")` |
| `LaBrRun.histFile` (nts .root): `DigiLaBr/ttree` | per sample: eventId, time, ADC, chargeCollected, chargeDecayed | `src/plotLaBrWaveform.C`, `src/compareLaBrSpectra.C` |
| `TruthLaBr/ttree`, `TruthLaBr/decays` | truth energy per event; per decay (energy, source, time) | `src/compareLaBrSpectra.C`, `src/plotIntrinsicTruthSpectrum.C` |
| `TagPileup/pulses` | truth pulses on the absolute timeline with pileup tags | `plotLaBrBinary.py --truth` (after exporting to CSV, see `--help-truth`) |
| `LaBrRun.artFile` (.art) | every product | re-run truth/pileup without Geant4: `fcl/LaBrTruthPulsesFromArt.fcl`; `src/dumpLaBrTruthPulses.C` |

Mapping truth to waveform samples, per event:
sample = (eventId - 1) x 628 + floor(t_local / 2.7), where t_local = t_abs - (eventId - 1) x 1695.
The waveform peak is about 18 samples later.

---

## 9. Parameters that must be kept in sync

| If you change | Also change | Why |
|---|---|---|
| `LaBrPipeline.DigiLaBr.peakKeVPerADC` | `LaBrPipeline.TagPileup.peakADCPerKeV` (= 1/value); `PEAK_KEV_PER_ADC` in `plotLaBrBinary.py` | pileup windows and plot energy axes |
| `LaBrPipeline.DigiLaBr.pulseTau2` | `LaBrPipeline.TagPileup.pulseTau2` | pileup window length |
| `LaBrPipeline.DigiLaBr.fADC` | `FADC_MHZ`, `SAMPLES_PER_EVENT` in `plotLaBrBinary.py`; the pulse plot in `plotLaBrWaveform.C` (370.37 hard-coded) | sample count and spacing |
| `alphaQuenchingFactor` / `Slope` on DigiLaBr | the same keys on `LaBrPipeline.TruthPulses` and `.TruthLaBr` | truth energies must match the digitized ones |
| `lightYield`, `quantumEfficiency` | `plotLaBrWaveform.C` / `compareLaBrSpectra.C` pe -> keV constant (63000 x 0.35) | charge-based energy in those macros |
| `LaBrRun.nEvents` | nothing (it sets `TagPileup.nMerge` itself), but don't use `mu2e -n` for full runs | one pileup batch |
| any activity scale | every source's scale, if you compare rates between sources | the ratios |
| crystal geometry | everything in 7.3 | selection, decay placement |
| the micropulse length | not currently possible: the digitizer hard-codes 1695 ns | |

---

## 10. Validation status and known limitations

**What has been validated:**
- **Gun mode:** bit-identical to the old `PhotonGunToLaBrWaveform.fcl` (same photons, same
  seeds).
- **Poisson + La-138:** statistically consistent with the old
  `PoissonPhotonGunPlusLa138ToLaBrWaveform.fcl` (`compare.md`).
- **`spectrum`/`lines` + Poisson + disk beam** has run: the SSC FOILOFF 20 kHz job. Its
  digitized per-pulse spectrum matches truth, and truth pulses land on the waveform pulses.

**What has not been run yet:**
- `dataset` mode;
- `STMParticleDump`;
- the `histogram` and `rootHist` spectrum types;
- the cone option.

Treat these as untested until you've checked their outputs.

**Limitations:**
- **Digitizer RNG not seeded:** the same photostatistics sequence in every job
  (`STM_Simulation_Overview.md` section 11).
- **`NoiseSD` is not a working PMT noise model.** The value in mV is converted with the HPGe
  preamp feedback capacitance into "charge carriers" and then added to a signal in
  photoelectrons. The prolog's `PreamplifierNoiseSD : 0.32` would give sigma of about 3900 ADC.
  Keep `NoiseSD : 0` until a LaBr noise model exists.
- **Crystal geometry is hard-coded** in three C++ modules (Section 7), and there is the 101.6 mm
  z offset of unknown origin.
- **Micropulse fixed at 1695 ns**, with no beam spill/supercycle structure: every event is an
  in-spill micropulse.
- **Source classification is coarse:** everything from STMParticleSource is "external" (source
  0). `primary` and `background` cannot be separated in truth.
- **Timing:** one timing model per source. No time-within-micropulse profile (e.g. prompt X-rays
  vs the delayed 1809) other than `fixed` / `uniform` / `dataset`.
- **Dataset `transform`** has no special handling of the HPGe-absorber shadow (unlike the
  production shift), and it focuses angles by the radial scale.
- **Resolution is an input, not a prediction.** F = 1.27 was tuned at 662 keV to a catalogue
  value. At 1809 keV the model is about 15% more optimistic than the ELBE measurement.
- **Not modelled:** rate-dependent PMT gain drift and baseline noise.

---

## 11. Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `art ... exit status 89, No configuration file given` | You typed `mu2e build`; it's `muse build -j8` |
| `unknown module type STMParticleSource` / old behaviour after a rebuild | The build failed to link, so the old `.so` is still in use. `find build -name "*STMParticleSource*.so"` and read the build log |
| error about `process_name` | underscores are not allowed; `makeLinesFcl.C` strips them for you |
| FHiCL parse error on a line that looks fine | a non-ASCII character (em dash, smart quote) |
| `mode "gun" needs E > 0` / `needs a spectrum table` / `needs a dataset table` / `modes gun/spectrum need a beam table` | a table was replaced by one missing that key. Check the override order (3.1) |
| `emission.rateHz needs poisson : true` | you set a rate with `single` emission |
| `no dataset particles pass the selection` | wrong tree name or branch names, filters too tight, or the transform's `inputRadius` removes everything. The log line gives entries read, kept and filtered |
| `cannot read TTree ... from ...` | the file path is wrong. Relative paths are relative to the run directory, else `MU2E_SEARCH_PATH` |
| digitizer `rejected: out-of-bounds` large | the crystal constants don't match the Geant4 geometry (Section 7) |
| `rejected: step.time > buf` large | particles arrive more than 2 micropulses after the event start: raise `microspillBufferLengthCount`, or use `timing : "dataset"` |
| empty `TagPileup/pulses` | `TruthPulses`/`TagPileup` not in `truthProducers`, or an old build (before the end-of-job flush) |
| truth and waveform don't line up in your own code | you used a global time -> sample mapping; map per event (Section 8) |
| spectrum lines look 10 keV wide / offset | you used `rootHist` on a histogram with 10 keV bins. Use `lines` (or `makeLinesFcl.C`) |
| `LaBrDigi WARNING: ADC value ...` | a positive sample above 5000: shouldn't happen with the current signal sign. Report it |
