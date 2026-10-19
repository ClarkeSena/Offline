# Mu2e STM Simulation — Full Technical Overview

**Scope:** everything in `Offline/STMMC/` — what the simulation does, how each stage works,
the physics behind the two digitizers, the truth/pileup machinery, where the work currently
stands, and the traps that have cost time before.

**Written:** 2026-09-21, reconstructed from the code as it stands on disk.

**Read this first.** The other documents in this directory are narrower and older:

| Document | Covers |
|---|---|
| `La138_STATUS.md` | La-138 intrinsic-activity work — status and file map |
| `CalibrationSource.md` | Known-activity source runs; validating the La-138 rate |
| `La138_MethodC_implementation.md` | Full derivation of the La-138 decay-product generator |
| `PhotonGunPlusLa138_TruthSpectrum.md` | Truth energy spectrum split by source |
| `FullSim_HPGe_LaBr.md` | Combined two-detector simulation |
| `STM_3D_Visualization.md` | Geometry viewing via ROOT TGeo + GDML |
| `../../summary.md` (repo root) | Handoff note from the ELBE detector-response update |

---

## 1. What the simulation is for

The **Stopping Target Monitor (STM)** measures the rate of muons stopping in the Mu2e
stopping target, by counting characteristic X-rays and gammas emitted from it. That rate is
the denominator of the conversion-rate measurement, so the STM is a normalisation detector.

Two detectors sit at the downstream end of the beamline:

- **HPGe** — high-purity germanium, a semiconductor diode. Superb energy resolution
  (sub-keV), slow signal formation.
- **LaBr₃:Ce+Sr** — a scintillator crystal read out by a Hamamatsu R10233-100 PMT. Much
  faster, worse resolution, tolerates higher rates.

This package simulates the **electronics response** of both: it takes Geant4 truth energy
deposits and produces realistic digitized ADC waveforms that look like what the real DAQ will
record. The point is to develop and validate the offline reconstruction — pulse finding,
pileup handling, energy calibration — against a dataset where the truth is known exactly.

### The central design principle

Geant4's job in this pipeline is **truth energy deposits only**: energy, position, time.

Mu2eG4 has **no optical physics and no `MaterialPropertiesTable` for any Mu2e detector** —
not the STM, not even the CRV. This is deliberate, not an oversight. Simulating individual
optical photons is expensive and requires measured optical constants (reflectivity, bulk
absorption length, refractive indices, surface finish, wrap) that generally do not exist for
the as-built detector.

Therefore **all** light production, PMT response, pulse shape and electronics behaviour is
**parametric, inside the digitizer module**. The CRV takes the same approach with a
parametrised photon-yield lookup.

One asymmetry between the two detectors is worth internalising, because it explains a lot of
the code:

- For **LaBr**, scintillation light is essentially position-independent, so the deposit
  position is used *only* as an acceptance cut. The pulse shape is a fixed template.
- For **HPGe**, signal formation genuinely *is* position-dependent — the induced current
  depends on where the charge carriers start and how far they drift. So the HPGe digitizer
  derives the pulse rise from drift physics, and the deposit position feeds directly into
  the waveform.

That asymmetry is physically justified, not an inconsistency.

---

## 2. Framework and time structure

Everything runs inside **art**, the Mu2e/FNAL event-processing framework, configured through
FHiCL files in `Offline/STMMC/fcl/`.

The single most important convention:

> **One art event = one micropulse = 1695 ns.**

The Mu2e beam is delivered in micropulses 1695 ns apart. A long run is therefore a
*concatenated stream* of micropulses, and several modules maintain running state across
events to stitch them into one continuous timeline:

- The digitizers carry **charge carry-over** buffers, so a pulse late in one micropulse
  spills correctly into the next.
- `LaBrTruthPulses` maintains a running **time origin** that advances by 1695 ns every
  event, putting all truth pulses on one absolute axis.
- `ConcatenateDigitizedWaveforms` and `TagPileupPulses` both buffer `nMerge` events and
  emit on the batch boundary.

At the LaBr sampling rate of 370.370370 MSPS (tADC = 2.7 ns), one micropulse is

```
nADCs = (int)(1695.0 / 2.7) + 1 = 628 samples
```

**Subtle but important:** 628 samples × 2.7 ns = 1695.6 ns, but the truth timeline advances
exactly 1695.0 ns per event. The two drift apart by 0.6 ns per event — about 2 samples over a
10,000-event run. Any mapping from truth time to a sample index in the concatenated stream
**must be done per event**, not with a single global `time × fADC`. `dumpLaBrTruthPulses.C`
does this correctly; anything new must too.

---

## 3. The pipeline

```
 GENERATION            GEANT4                DIGITIZATION           OUTPUT
 ──────────            ──────                ────────────           ──────
 PhotonGun         →   g4run              →  LaBrWaveforms...    →  STMWaveformDigi
 PoissonPhotonGun      QGSP_BERT_EMZ         HPGeWaveforms...       nts.*.root TTree
 La138Source           SDConfig: STMDet      (parametric            labr_*_ADC.bin
 GenParticleMerger     → StepPointMC          response)             (raw int16 stream)

 TRUTH (parallel branch)
 ───────────────────────
 g4run:STMDet → LaBrTruthSpectrum                       → TruthLaBr/ttree
 g4run:STMDet → LaBrTruthPulses → TagPileupPulses       → STMMCPulse + pulses/ttree
```

### 3.1 Generation

| Module | Purpose |
|---|---|
| `PhotonGun` (EventGenerator) | Fixed-energy photons from a fixed point. Simple tests. |
| `PoissonPhotonGun_module.cc` | Realistic Poisson pulse train with a configurable spectrum. |
| `La138Source_module.cc` | Intrinsic LaBr radioactivity (see §5). |
| `GenParticleMerger_module.cc` | Concatenates several `GenParticleCollection`s into one for `g4run`. |

**`PoissonPhotonGun`** is the realistic source. Per event it fires `Poisson(meanPhotonsPerEvent)`
photons, each at an independent uniform time in `[0, maxTime)`. Concatenating micropulses then
gives a genuine homogeneous Poisson process of rate `meanPhotonsPerEvent / maxTime`:

```
0.0339 photons / 1695 ns = 20 kHz
```

At 20 kHz roughly 96.6% of events are empty, so realistic runs need many events — 590,000
events is one second of beam.

Each photon's energy is drawn from a configurable spectrum:

- a single fixed line (`E`), and/or
- several discrete lines (`lineEnergies` with relative `lineIntensities`), and/or
- an exponential continuum (`background: { enabled, weight, E0, EMin, EMax }`), pdf
  `~ exp(-E/E0)` truncated to `[EMin, EMax]`, sampled by inverse CDF.

A photon is assigned to a component in proportion to its weight. The intensities set only the
*mix*; the total rate stays `meanPhotonsPerEvent`. There is a guard against
`pt > min(E)` (which would make `pz = sqrt(E² - pt²)` imaginary).

When merging generators, remember to point Geant4 at the merged product:

```
physics.producers.g4run.inputs.primaryTag: "genCombined"
```

### 3.2 Geant4

`g4run` with `physicsListName: "QGSP_BERT_EMZ"` and `SDConfig.enableSD: [STMDet]`.

All the real spectral structure originates here, from honest transport: the full-energy
photopeak, the Compton edge, single- and double-escape peaks, backscatter, and coincidence
sums. For a 1809 keV gamma you should expect the Compton edge at ~1585 keV, single escape at
~1298 keV and double escape at ~787 keV. **A continuum below the photopeak is correct
physics, not a bug.**

Output is a `StepPointMCCollection` tagged `g4run:STMDet`.

### 3.3 Digitization

Both digitizers read the **same** `g4run:STMDet` collection and self-select their own crystal.
There is no separate per-detector step collection.

### 3.4 Output

- `STMWaveformDigi` — the art data product, one per micropulse, holding a trigger time offset
  and the `int16_t` ADC vector.
- `nts.*.root` — TTrees (`DigiLaBr/ttree`, `DigiHPGe/ttree`, `TruthLaBr/ttree`,
  `TagPileup/pulses`) for plotting.
- `BinaryWaveformWriter_module.cc` — writes the ADC samples as a **headerless little-endian
  int16 stream**, waveforms concatenated in event order. This is the raw-stream format the
  offline matched-filter pulse finder consumes. Read it with
  `np.fromfile(fname, dtype="<i2")`.

---

## 4. The two digitizers

| | HPGe | LaBr₃:Ce+Sr |
|---|---|---|
| Module | `HPGeWaveformsFromStepPointMCs_module.cc` | `LaBrWaveformsFromStepPointMCsTest_module.cc` |
| Detector type | Semiconductor diode | Scintillator + PMT |
| Signal origin | Charge-carrier drift (Shockley-Ramo) | Scintillation light → photoelectrons |
| Readout model | Charge-sensitive preamp, IIR filter | Direct PMT anode, no integration |
| Statistics | Fano factor 0.1, 2.96 eV/pair | Poisson N_pe + PMT excess noise F = 1.27 |
| Rise time | 220–513 ns, position-dependent (drift) | 23 ns 10–90%, fixed (fitted SPE template; peak at 49.95 ns) |
| Decay | ~50 µs preamp RC | ~28 ns tail constant |
| fADC | 320 MSPS (3.125 ns, 543 samples) | 370.370370 MSPS (2.7 ns, 628 samples) |
| Crystal | Coaxial, bored, rotated −45° about y | Solid cylinder, identity rotation |
| Product instance | `produces<...>("HPGe")` | `produces<...>("LaBr")` |

Both produce **negative-going** pulses and clip at `ADCMax = -32767` (`-2^15 + 1`), since the
output is `int16_t`.

### 4.1 LaBr signal chain, step by step

**(a) Crystal selection.** LaBr sits on the +x side of the STM. Steps with `x < -3925` are
rejected as HPGe-side. Then the step is transformed to crystal-local coordinates (identity
rotation, so just subtract the crystal centre and shift the front face to z = 0) and required
to lie inside the solid cylinder, `R ≤ 38.2`, `0 ≤ z ≤ 76.3` (radius/length plus a 0.1 mm
tolerance). There is also a time cut: `step.time() > microspillBufferLengthCount × 1695 ns`
is dropped.

**(b) Energy → photoelectrons.**

```
N_pe_mean = E_dep × lightYield × quantumEfficiency
          = E_dep[MeV] × 63000 × 0.35
          = 22050 pe per MeV
N_pe      = Poisson(N_pe_mean)
```

**(c) PMT excess noise.** Dynode-chain multiplication broadens the measured charge beyond
pure photostatistics. The digitizer adds Gaussian excess variance `(F-1)·N_pe` on top of the
Poisson count, giving total relative variance `F/N_pe`:

```
signal_pe = N_pe + Gauss(0,1) × sqrt((F - 1) × N_pe)
```

This is the standard photostatistics-with-excess-noise result, `σ/E = sqrt(F/N_pe)`.
**F = 1.27 is not a free guess** — it is calibrated so that with QE = 0.35 the model
reproduces the **2.2% FWHM at 662 keV** measured in the ELBE test-beam qualification. Poisson
alone would give 1.95%. The assumption baked in is that all the sub-Poisson broadening is PMT
excess noise. At 1809 keV this gives roughly 24 keV FWHM.

**(d) Pulse shape.** The single-photoelectron response is

```
shape(t) = (1 - e^{-t/τ₁})^p · e^{-t/τ₂}
τ₁ = 22.3543 ns    τ₂ = 28.0007 ns    p = 6.6596
```

fitted to real LaBr₃ data. The same function is the matched-filter template used by the
pulse finder on real data, which is exactly why it lives here in this form.

Because *p* is non-integer the shape has **no closed-form integral**, so the constructor
tabulates its normalised CDF by the trapezoidal rule on a 0.05 ns grid out to 40·τ₂ (the tail
beyond that is of order `e^{-40}`, negligible) and interpolates linearly. Each ADC bin then
receives

```
fraction = pulseCDF(dt + tADC) - pulseCDF(dt)
```

so the bins sum to 1 and **charge is exactly conserved**.

The peak of the shape has a closed form, `t_peak = -τ₁·ln(τ₁/(p·τ₂ + τ₁)) = 49.95 ns`
≈ 18.5 samples. So the waveform maximum arrives ~50 ns *after* the truth deposit time. This
matters constantly when comparing truth to waveform.

The electronics rise is the PMT anode + amplifier single-electron response. Geant4 does not
model it, so it correctly belongs in the digitizer as a convolution.

**(e) Readout.** `decayCharge()` is a **pure passthrough** — a copy. The PMT anode voltage
tracks the photoelectron arrival rate directly; there is no charge-sensitive preamplifier to
integrate. This is the single biggest structural difference from HPGe, and the reason
`risingEdgeDecayConstant` (50 ns) is configured but currently unused in the LaBr signal path.

**(f) Noise.** `addNoise()` adds Gaussian noise of width `NoiseSD`, converted from mV to
charge-carrier units. **Currently `NoiseSD = 0` in every FCL** — see open items.

**(g) Digitization.** `ADC = charge × chargeToADC`, rounded, clipped at `ADCMax`.

### 4.2 Gain calibration — two modes

This trips people up, so be explicit about which mode a job is in. The `beginJob()` printout
with `verbosityLevel: 1` states it (`calibration mode  PEAK` or `INTEGRAL`).

**INTEGRAL mode** (`EnergyPerADCBin`, keV/bin, default 0.5) makes the **pulse integral** the
energy estimator:

```
chargeToADC = (1e6 / pePerMeV) / (EnergyPerADCBin × 1e3)
```

**PEAK mode** (`peakKeVPerADC`, keV/ADC — prolog sets **0.57** from source data) makes the
**peak height** the estimator, `E = peakKeVPerADC × ADC_peak`:

```
chargeToADC = (1e3 / pePerMeV) / (peakKeVPerADC × peakFraction)
```

`peakKeVPerADC > 0` **overrides** `EnergyPerADCBin`.

`peakFraction` is the share of a single deposit's signal that lands in its tallest bin. The
digitizer **computes it from its own CDF table** at construction, by scanning all bins for
the maximum per-bin fraction. That means the peak calibration stays valid if `fADC` or the
pulse-shape parameters change — a genuinely nice property. Values: **0.04114** at 370.37 MHz,
**0.05065** at 300 MHz.

At `peakKeVPerADC = 0.57`: 1809 keV → 3174 ADC peak, 662 keV → 1161 ADC. The `ADCMax` floor
clips above roughly 13 MeV.

**Known caveat:** deposits start on bin boundaries (`tStart` is an integer division of
`step.time()` by `tADC`), so there is **no sub-sample jitter**. Simulated peak-height
resolution is therefore cleaner than real source data. If you ever need to match real peak
resolution exactly, this is the first thing to fix.

### 4.3 HPGe signal chain, in brief

Physically different enough to be worth stating:

**Geometry.** Crystal centre `(-3986.30, 0, 40612.70)` mm, length 78.5 mm, radius 36.05 mm,
with a coaxial bore hole (radius 5.25 mm, length 64.7 mm plus a hemispherical end). Steps are
transformed by `rotateY(-45°)` so the crystal axis points along +z. The module prefers
GeometryService values and falls back to the hard-coded constants; it prints the difference if
they disagree.

> The `-45°` sign matters. The original code applied `rotateY(+45°)` — the forward rotation
> instead of the inverse — which was wrong.

**Charge creation.** `N_eh = E_dep / 2.96 eV`, with Fano fluctuation
`N += Gauss(0,1)·sqrt(F·N)`, F = 0.1 (Ge at 77 K). Per-step fluctuations add in quadrature,
so the total variance over all steps is `F·N_total`, as required.

**Charge collection (Shockley-Ramo).** Electrons and holes drift at approximately 0.08 and
0.06 mm/ns respectively. For a coaxial geometry the weighting potential is logarithmic, so
the induced charge as a function of time is

```
Q(t) = A · ln[ (1 + B_h·t) / (1 - B_e·t) ],   A = N_eh / ln(R₂/R₁)
```

with separate branches for the phases after one carrier species has finished drifting.
**Rise time is therefore position-dependent** — a deposit near the bore sees a different rise
from one near the outer surface. This is the physical reason HPGe needs the position and LaBr
does not.

**Preamplifier.** A charge-sensitive preamp modelled as a one-pole IIR filter with a ~50 µs
decay:

```
chargeDecayed[t] = chargeDecayed[t-1] · exp(-tADC/τ) + chargeCollected[t]
```

carrying `lastEventEndDecayedCharge` across event boundaries, with a periodic reset at
`resetEventNumber` (31858) mimicking the real preamp reset.

### 4.4 The detector-selection x cuts

HPGe rejects `x > -3904`; LaBr rejects `x < -3925`. Those two windows **overlap** in the band
`-3925 < x < -3904`, so in principle a step there would be accepted by both.

In practice the band is empty: after the crystal-cylinder bounds check, HPGe only accepts
x ∈ [-4022, -3950] and LaBr only x ∈ [-3901, -3825]. The geometry cut disambiguates, so this
is harmless — but it is worth knowing it is the *cylinder* cut, not the x cut, doing the real
work.

---

## 5. La-138 intrinsic radioactivity

Real LaBr₃ spectra show peaks the simulation originally lacked, because the simulation only
fired external particles. Lanthanum has a naturally radioactive isotope, **¹³⁸La** (0.09%
abundance), and the crystal is therefore its own radioactive source:

- **EC branch (66.4%)**: 1436 keV gamma, plus (71% of the time) a Ba K X-ray at 32.2/36.4
  keV. Because the source is *inside* the detector, these sum → the characteristic
  **1468 keV sum peak** with a **1436 keV escape peak** alongside.
- **β⁻ branch (33.6%)**: 789 keV gamma plus a beta electron with a 263 keV endpoint. Again,
  internal, so the two sum → a **789–1052 keV continuum**, *not* a sharp 789 keV line. That
  is the correct answer, and it looks wrong at first glance.

The activity is **computed from the material**, not hard-coded:

```
N(La-138) = mass / molarMass × N_A × laPerMolecule × la138Abundance
activity  = ln2 / halfLife × N(La-138)
```

giving ~540 Bq (1.57 Bq/cm³) for the 348 cm³ crystal — about 9 × 10⁻⁴ decays per 1695 ns
event. You need ~1M events to see the peak at the physical rate, so `activityScale` (default
1.0) multiplies the rate for statistics studies; the combined FCL uses 1000.

### Why this needed "Method C" — the hard-won lesson

Two earlier approaches both failed, for reasons worth remembering:

1. **Let Geant4 decay the crystal.** It will not. `G4RadioactiveDecay` acts on ion
   **tracks**, never on the atoms making up a material. The material only determines how many
   decays there are and where. You must inject ions as primaries.
2. **Inject ¹³⁸La ions as primaries.** Still zero decays. Analogue `G4RadioactiveDecay`
   samples a decay time from the half-life — 1.02 × 10¹¹ years. The in-window decay
   probability is ~10⁻²⁵, so the sampled time never arrives inside a microsecond event. (A
   biased-decay physics list, `Mu2eG4BiasedRDPhysics` via `radiationVRmode`, would force
   decays in a time window with weights — a viable but heavier alternative.)

   Along the way: the G4 ion PDG code is `1000571380`, from `1e9 + 10000·Z + 10·A` with
   Z = 57, A = 138. An earlier `1000570138` was simply wrong.

3. **Method C (adopted).** Emit the **decay products directly** as `GenParticle`s at a
   uniform vertex inside the crystal cylinder, sharing position and time. Geant4 then
   transports them normally with `turnOnRadioactiveDecay: false`, so all the photopeak,
   Compton, escape and sum structure emerges from honest transport, and events are
   unit-weight.

The generator was standalone-validated: branching 0.664/0.336, X-ray-per-EC 0.71, the 1468
sum and 1436 escape peaks, and the 789–1052 continuum all come out correctly.

---

## 6. Truth machinery

### 6.1 Truth energy spectrum — `LaBrTruthSpectrum_module.cc`

An `EDAnalyzer` that builds the Geant4 truth energy spectrum from `g4run:STMDet`, split into
photon-gun signal vs La-138 intrinsic activity. Per-event sums go into `TruthLaBr/ttree` as
`edepTotal`, `edepGun`, `edepLa138` (keV).

Classification traces each step back to its true primary and thresholds on the primary's
**start z** (`photonGunMaxZ`, default 40700 mm): the gun fires from the aperture at z = 40404,
La-138 originates inside the crystal at z ≈ 40891.

Truth and digitized spectra differ **only** by detector resolution — they share the same
22050 pe/MeV calibration — so overlaying them bin-for-bin isolates exactly what the digitizer
smearing does.

> ### The `originParticle()` trap
>
> `SimParticle::originParticle()` does **not** walk to the primary. It only recurses for
> `selfParent()` — that is, when `hasParent()` and `creationCode == mu2ePrimary`. For an
> ordinary secondary (a Compton or photoelectric electron) it returns **the secondary
> itself**, whose vertex is inside the crystal. Every step then classifies as La-138.
>
> You must walk `parent()` manually until you reach the particle with no parent:
>
> ```cpp
> const SimParticle* primary = step.simParticle().get();
> while (primary && primary->hasParent()) {
>   const SimParticle* up = primary->parent().get();
>   if (!up) break;
>   primary = up;
> }
> ```
>
> Both `LaBrTruthSpectrum` and `LaBrTruthPulses` do this, identically and deliberately.

### 6.2 Pileup truth — the two-stage chain

This mirrors the digitize → concatenate architecture, and produces a labelled dataset for
evaluating the offline pulse finder.

**Data product:** `mu2e::STMMCPulse` / `STMMCPulseCollection` in
`MCDataProducts/inc/STMMCPulse.hh`, registered in that package's `classes.h` and
`classes_def.xml`.

A **truth pulse** is: the total ionizing energy **one primary** deposits in the crystal,
placed at its **earliest** deposit time (the leading edge the pulse finder triggers on), on an
**absolute** time axis.

Fields: `time`, `energy`, `primaryId`, `source` (`unknown`/`photonGun`/`La138`), `eventId`,
plus the tags filled in stage 2 — `isPileup`, `nOverlap`, `dtNearest`, `window`.

**Stage 1 — `LaBrTruthPulses_module.cc`.** Per micropulse, groups `g4run:STMDet` steps by
traced primary and emits one `STMMCPulse` each. Crystal selection and primary tracing match
`LaBrTruthSpectrum` exactly, so the pulses cover the same volume the digitizer sees.
`tAbs = _timeOrigin + tMin`, where `_timeOrigin` advances by `micropulseTime` **every event,
whether or not it had a deposit** — so the input file must be read whole and in order.
`isPileup` is left false here.

**Stage 2 — `TagPileupPulses_module.cc`.** Buffers `nMerge` micropulses onto one sorted
timeline and tags overlaps. Working on the merged sorted list makes cross-micropulse pileup
automatic: the last pulse of one event and the first of the next are simply adjacent entries.

The window is energy-dependent and analytic. A pulse of energy E has peak amplitude
`peakADCPerKeV × E`, and its tail decays as `exp(-dt/τ₂)`, so it stays above the resolving
threshold for

```
w_after(E) = τ₂ · ln( peakADCPerKeV × E / thresholdADC )     [floored at minWindow]
```

Two pulses overlap if `dt < w_after(earlier)` (the later one sits on the earlier one's
still-sloping tail) **or** `dt < wBefore` (a fixed rise-side window, 20 ns, where the earlier
pulse is caught on the later one's rising edge).

Because `w_after` grows only **logarithmically** with energy, the reach past a micropulse
boundary is at most a few hundred nanoseconds even for very large pulses — so buffering a
batch captures essentially all cross-boundary pileup. Cross-microspill energy dependence is
handled by carrying the `(tAbs, E)` tuple in the buffer and **recomputing** `w_after`, not by
carrying charge forward.

> **Set `nMerge` = the total event count** so the whole run is one batch. Otherwise pileup is
> lost at every batch boundary — the same limitation waveform concatenation already has.

### 6.3 Second pass without Geant4 — `LaBrTruthPulsesFromArt.fcl`

`LaBrTruthPulses` needs only the `StepPointMC`s and the `SimParticle` parent chain, both of
which a plain `RootOutput` already stored. So any job that wrote `g4run:STMDet` can be
reprocessed to rebuild and re-tag truth pulses **without re-running Geant4**. This is what
made iterating on the pileup-window parameters cheap.

Two requirements: `process_name` must **differ** from the input file's process name, and
`nMerge` must equal the input file's event count.

### 6.4 Truth ↔ waveform mapping — `dumpLaBrTruthPulses.C`

Dumps truth pulses to flat CSV **and locates each one on the concatenated ADC `.bin`**. It
reads the art file through ROOT's emulated streamers, so it needs no Offline dictionaries and
runs under a bare ROOT anywhere — not just on the gpvm.

The index conventions, all validated (96.3% of isolated >100 keV pulses land exactly on the
predicted `peakSample`, 100% within one sample):

- Event N's 628 samples occupy `.bin` indices `[(N-1)·628, N·628)`; sample j of event N is at
  absolute time `(N-1)·1695 + j·2.7` ns. **Per-event**, because of the 0.6 ns/event drift.
- `STMMCPulse::time` is the **leading edge**, not the peak. The maximum arrives
  `t_peak ≈ 49.95 ns` later.
- ADC bins integrate charge over the bin, so a feature at time t lands in the bin whose
  *centre* is nearest t: `floor(t/tADC - 0.5)`.
- Pulses are negative-going; `peakADC` is the magnitude.

Output columns: `pulseIndex, eventId, primaryId, source, sourceLabel, edgeTime_ns,
localTime_ns, edgeSample, peakTime_ns, peakSample, inBin, carryOver, energy_keV, peakADC,
isPileup, nOverlap, dtNearest_ns, window_ns`.

---

## 7. Key parameters — `fcl/prolog.fcl`

### `LaBrDigitization`

| Key | Value | Meaning |
|---|---|---|
| `EnergyPerADCBin` | 0.5 keV/bin | Integral-mode calibration |
| `peakKeVPerADC` | 0.57 keV/ADC | **Peak-mode calibration; overrides the above** |
| `pulseTau1` | 22.3543 ns | SPE rise constant τ₁ |
| `pulseTau2` | 28.0007 ns | SPE decay constant τ₂ |
| `pulseRiseSharpening` | 6.6596 | Rise sharpening exponent p |
| `lightYield` | 63000 photons/MeV | LaBr₃:Ce+Sr scintillation yield |
| `quantumEfficiency` | 0.35 | R10233-100 super-bialkali QE |
| `excessNoiseFactor` | 1.27 | PMT dynode excess noise F |
| `DecayConstant` | 0.05 µs | Preamp decay — **configured but unused** (LaBr has no preamp stage) |
| `PreamplifierNoiseSD` | 0.32 mV | Not currently wired in; jobs set `NoiseSD: 0.0` |

### `LaBrPileup`

| Key | Value | Meaning |
|---|---|---|
| `micropulseTime` | 1695.0 ns | Must match the digitizer |
| `photonGunMaxZ` | 40700.0 mm | Primary start-z split: below = gun, above = La-138 |
| `pulseTau2` | 28.0007 ns | Tail constant, mirrors `LaBrDigitization` |
| `peakADCPerKeV` | 1.7543859649 | = 1 / 0.57 |
| `thresholdADC` | 5.0 | Resolving threshold; should match the pulse finder |
| `wBefore` | 20.0 ns | Fixed rise-side window |
| `minWindow` | 20.0 ns | Floor on the forward window |

### `HPGeDigitization`

| Key | Value |
|---|---|
| `EnergyPerADCBin` | 0.5 keV/bin |
| `FanoFactor` | 0.1 (Ge at 77 K) |
| `DecayConstant` | 50 µs |
| `resetEventNumber` | 31858 |
| `concatenation.nMerge` | 31858 |
| `concatenation.STMWaveformDigisTag` | `DigiHPGe:HPGe` |

### `STMDAQParameters.samplingFrequencies`

HPGe 320 MSPS; LaBr 370.370370 MSPS.

### Other stanzas

`ComponentPositions` (beam-order STM component positions), `MixSTMEvents` (Poisson mixing
rates for production), `STMMCAnalysis` (ZS/MWD reco tags), `SimplifyStage1Data` /
`SimplifyStage2Data`, `ResamplingProducer` / `ResamplingFilter`, `ShiftVD101Steps`,
`Efficiency`.

---

## 8. Geometry

Mu2e global coordinates, mm.

| | Value | Notes |
|---|---|---|
| HPGe crystal centre | (-3986.30, 0, 40612.70) | Coaxial, bored, `rotateY(-45°)`; L = 78.5, R = 36.05, bore R = 5.25 / L = 64.7 |
| LaBr crystal centre | (-3863.36, 0, 40891.50) | Solid cylinder, identity rotation; L = 76.2, R = 38.1 |
| HPGe aperture | (-3944.6, 0, 40404) | Where the photon gun fires from |
| LaBr aperture x | -3863.4 | |
| Absorber | (-3944.6, -24.9, 39900) | 50 × 50 × 390 mm; the y offset is an open question in prolog |
| VD101 | (-3904.0, 0, 40301.0) | |

> ### The 101.6 mm z-offset
>
> `STMGeometryDump_module.cc` reports crystal z positions **101.6 mm higher** than the actual
> G4 step positions, for **both** detectors (HPGe 40714.3 → 40612.7; LaBr 40993.1 → 40891.5).
> Both digitizers use the empirically corrected values. **The root cause is still unknown.**
> 101.6 mm = 4 inches exactly, which smells like a duplicated offset somewhere in the
> geometry construction, but that has not been chased down.

> ### Documentation inconsistency to be aware of
>
> `CLAUDE.md` lists the HPGe crystal centre as **(-3942.15, 0, 40612.70)**. The module uses
> **(-3986.30, 0, 40612.70)**. The -3942.15 value is `ComponentPositions.HPGeCrystalEndcapCentre`
> — the endcap, not the crystal centre. Trust the module.

**Viewing the geometry:** this Muse build (`al9-prof-e29-p092`) is compiled **without Geant4
visualisation**. In `Mu2eG4_module.cc` the entire vis path — creating `G4VisManager` *and*
executing `visualization.initMacro` — sits behind `#if (defined G4VIS_USE_OPENGLX || ...)`, so
any `.mac` file is silently ignored: no window, no `/vis/list`, no file export. The per-event
"Enter a character to go to the next event" prompt still appears but does nothing. Do not
chase OGLSX/OGLSQt in this build.

Instead, export GDML (`g4run.debug.writeGDML: true` + `GDMLFileName`) and render it in ROOT's
TGeo OpenGL viewer, which negotiates GL differently and works over XQuartz.
`mu2e_common.gdml` is already in the working directory; `src/viewSTMClean.C` is the clean-view
macro. Logical volume names: `TungstenSSCPV` (SSC collimator), `fLaBrPV`, `fHPGePV`,
`stmDownstreamEnvelope` (STM mother).

---

## 9. File map

### Digitizers and generators

| File | Role |
|---|---|
| `src/LaBrWaveformsFromStepPointMCsTest_module.cc` | **The LaBr digitizer.** Main working module. |
| `src/HPGeWaveformsFromStepPointMCs_module.cc` | **The HPGe digitizer.** The promoted, fixed-geometry version. |
| `src/HPGeWaveformsFromStepPointMCs_module_test.cc` | Redundant — content was copied over the `_module.cc`. |
| `src/LaBrWaveformsFromStepPointMCs_module.cc` / `_module_test.cc` | Older LaBr attempts, superseded. |
| `src/PoissonPhotonGun_module.cc` | Poisson pulse train with configurable spectrum. |
| `src/La138Source_module.cc` | ¹³⁸La intrinsic activity (Method C). |
| `src/CalibrationSource_module.cc` | Sealed calibration source at a distance, absolute Bq rate. |
| `src/GenParticleMerger_module.cc` | Merges generator collections. |

### Truth and analysis

| File | Role |
|---|---|
| `src/LaBrTruthSpectrum_module.cc` | Truth energy spectrum, split by source. |
| `src/LaBrTruthPulses_module.cc` | Pileup truth stage 1 — one pulse per primary. |
| `src/TagPileupPulses_module.cc` | Pileup truth stage 2 — merge and tag overlaps. |
| `src/BinaryWaveformWriter_module.cc` | Raw int16 ADC stream for the pulse finder. |
| `src/STMGeometryDump_module.cc` | Dumps geometry positions (note the 101.6 mm offset). |
| `src/DetectorEfficiency_module.cc` | Photon interaction efficiency studies. |
| `src/ConcatenateDigitizedWaveforms_module.cc` | Concatenates `nMerge` micropulse waveforms. |

### Plot macros

| File | Produces |
|---|---|
| `src/plotLaBrWaveform.C` | Full waveform, individual pulses, energy spectrum (250 bins, 0–2000 keV) |
| `src/plotHPGeAndLaBrSpectra.C` | Overlaid HPGe + LaBr digitized spectra |
| `src/plotLa138TruthSpectrum.C` | Truth spectra — Combined / PhotonGun / La138 |
| `src/plotLaBrTruthPositions.C` | Truth deposit positions |
| `src/dumpLaBrTruthPulses.C` | Truth pulse CSV with `.bin` sample indices |
| `src/analyseCalibrationRates.C` | Absolute counts/s per source; source-to-La-138 ratio |
| `src/viewSTMClean.C` | ROOT TGeo geometry view |

The pre-ADC energy estimator used by the plot macros is the per-event sum of
`|chargeCollected|`, converted with `peToKeV = 1000 / (63000 × 0.35)` for LaBr and
`× 2.96/1000` for HPGe. It carries the detector resolution but **not** ADC quantization or
noise.

### Job configurations

| FCL | Pipeline |
|---|---|
| `PhotonGunToLaBrWaveform.fcl` | Simplest LaBr test: fixed 1.809 MeV gamma |
| `PhotonGunToHPGeWaveform.fcl` | Same for HPGe |
| `PoissonPhotonGunToLaBrWaveform.fcl` | **20 kHz Poisson train + full truth-pulse chain + binary writer** |
| `PoissonPhotonGunPlusLa138ToLaBrWaveform.fcl` | Poisson signal + intrinsic La-138 |
| `PhotonGunPlusLa138ToLaBrWaveform.fcl` | Fixed gun + La-138 |
| `La138OnlyToLaBrWaveform.fcl` | Intrinsic activity alone |
| `CalibrationSourceToLaBrWaveform.fcl` | Known source + La-138, for absolute rate validation |
| `LaBrTruthPulsesFromArt.fcl` | Second pass, no Geant4 |
| `FullSimToHPGeAndLaBrWaveforms.fcl` | Both detectors in one job |
| `HPGeReco.fcl` | Production reco: digitize → concatenate → ZS → MWD |
| `STMGeometryDump.fcl`, `VisSTMBeam.fcl` | Geometry dump / visualisation |

---

## 10. Current state of the work

The last stretch (July–August 2026) moved from *building* the digitizer to **producing a
labelled dataset for the offline pulse finder**. The chain in place:

1. `PoissonPhotonGun` generates a realistic 20 kHz multi-line source.
2. `BinaryWaveformWriter` dumps the concatenated ADC stream in the exact raw format the
   matched-filter pulse finder consumes.
3. `LaBrTruthPulses` → `TagPileupPulses` produce per-primary truth pulses with pileup labels
   on the same absolute timeline.
4. `dumpLaBrTruthPulses.C` maps each truth pulse onto a sample index in the `.bin`, validated
   to within one sample.
5. `LaBrTruthPulsesFromArt.fcl` re-runs steps 3–4 on an existing art file without Geant4, so
   the tagging parameters can be iterated cheaply.

### Artefacts on disk

| File | Contents |
|---|---|
| `PlotsFromGpvm/labr_sim_poisson_la138_ADC.bin` | 125.6 MB = **100,000 micropulses** (~170 ms of beam) of Poisson gun + La-138 at 20 kHz |
| `PlotsFromGpvm/sim.LaBrTruthPulses.art` | Output of the second-pass truth job (2026-08-07, newest) |
| `labr_la138_truth_pulses_v4.csv` | **6,622 truth pulses** — 637 photon-gun, 5,985 La-138 — of which **3,790 tagged pileup** |
| `labr_la138_truth_pulses{,_v2,_v3}.csv` | Earlier iterations: 4404/800, 6244/1978, 4610/1887 pulses/pileup |
| `15_07_labr_sim_20kHz_ADC.bin` + `15_07_labr_truth_pulses.csv` | The earlier 10k-event gun-only validation set (308 pulses) |

The v1 → v4 progression is the pileup tagging being tuned. The last thing in flight was
**iterating the truth-pulse and pileup labelling on the 100k-event dataset**, presumably so it
could be compared against what the pulse finder recovers from the `.bin`.

Note that at `activityScale = 1000`, La-138 pulses outnumber signal pulses roughly 9:1 in that
dataset — expected, but keep it in mind when quoting pileup fractions, since the fraction is
dominated by the scaled intrinsic rate rather than the physical one.

---

## 11. Open items

### Validation not yet done

- **Pileup fraction has not been checked against theory.** It should match the Poisson
  prediction `1 - exp(-R·w̄)` at 20 kHz. This is the obvious next check and it is cheap.
- **Overlay truth pulse times on the waveform, coloured by `isPileup`** — the visual
  sanity check that the labels mean what they should.
- **Compare against STMReco `peak_sep`** resolution, i.e. does the truth label agree with
  what the reconstruction can actually separate?
- **The combined HPGe + LaBr full sim has not been run end-to-end** on the gpvm; both 1809
  photopeaks and the LaBr intrinsic peaks should be verified to appear.

### Model refinements

- **Pileup window uses only the leading pulse's energy.** A two-energy version would let a
  large second pulse be resolved on a small first pulse's tail — currently it is tagged as
  pileup when in reality it would be separable.
- **Per-primary truth accumulation is per-event.** More than one La-138 decay in an event
  smears the 1436/789 lines together. Offered previously, not done.
- **`excessNoiseFactor` is anchored to the wrong number — the digitizer is ~15% too optimistic
  at 1809 keV.** The paper's two resolution figures are **not in the same units**: "2.2% FWHM
  at 662 keV" is FWHM and is a *cited Saint-Gobain catalogue spec* (ref [16]); "0.66% / 12 keV
  at 1808.7 keV" is **σ**, and is the paper's *own measurement*. F = 1.27 was solved for against
  the former. Converting to a common basis, the model gives σ = 10.2 keV at 1808.7 keV against
  a measured 11.9 keV. Three issues: (a) resolution depends only on `F/(lightYield·QE)`, so F is
  degenerate with the assumed 22050 pe/MeV and is not an independent measurement of dynode noise
  — though `F = δ/(δ-1)` does back-imply a plausible first-dynode δ ≈ 4.7; (b) the anchor is a
  catalogue figure at 662 keV rather than the measured value at the actual signal line; (c) a
  single multiplicative F forces `σ/E ∝ 1/√E` exactly, whereas real LaBr₃ has an
  energy-independent non-proportionality term. A two-term fit `(σ/E)² = a²/E + b²` to both paper
  values gives a constant floor b ≈ 0.43% σ. **Suggested fix:** keep F physical at ~1.27 and add
  an explicit `intrinsicResolution` term (~0.0043) in `depositCharge`, which reproduces both
  energies with every parameter keeping its meaning. Recalibrating F to 1.74 instead matches
  1808.7 keV but implies δ ≈ 2.4, turning F into a fudge factor. Note the paper explicitly
  endorses the *structure* (§5: "the measured energy resolution is fully accounted for by the
  photostatistics ... and the properties of the PMT") — only the anchor and the missing constant
  term are wrong.
- **Rate-dependent PMT gain (ELBE paper §6–8) is not modelled at all.** Photocathode
  charge-up causes gain *loss* (RC ~65 ms) while the last dynode causes gain *increase*
  (~25 ms equilibrium); the net is a few-percent drift over the 10–340 ms beam spill,
  corrected offline. This is a slow, beam-spill-level effect, not per-micropulse, so it needs
  a gain model keyed on a running-average singles rate across many microspills, plus the
  offline correction algorithm. Separate, larger feature.
- **Noise is switched off** (`NoiseSD = 0` everywhere). The HPGe charge-sensitive-preamp noise
  model is wrong for a PMT; LaBr needs a PMT-specific baseline noise model. This is distinct
  from `excessNoiseFactor`, which is gain-stage statistics, not baseline noise.
- **Alpha quenching and ²²⁷Ac contamination** are not modelled. Geant4 truth gives the full
  alpha energy and the digitizer applies one light yield to all particle types, so alphas land
  too high in energy. Needed for the alpha bumps seen in real data. ²²⁷Ac is also
  non-uniformly distributed, unlike ¹³⁸La.
- **Fermi function on the beta spectrum.** The allowed shape `N(E) ~ p·E·(E₀-E)²` is used with
  `F(Z=58, E)` omitted.
- **Sub-sample jitter.** Deposits start on integer bin boundaries, so simulated peak
  resolution is cleaner than real data.
- **Sampling rate worth confirming.** prolog uses 370.37 MHz for LaBr; the ELBE test beam used
  4 ns (250 MHz) 12-bit FADC. The Mu2e DAQ value has been left as-is.
- **`LaBrDigitization` prolog values** are a mix of ELBE-paper measurements and textbook
  placeholders. The three grounded in the paper are `quantumEfficiency`, the ~25 ns decay
  (now superseded by the fitted τ₁/τ₂/p) and `excessNoiseFactor`.

### Missing pieces

- **No HPGe truth-spectrum analyzer.** `LaBrTruthSpectrum`'s cylinder cut assumes identity
  rotation; HPGe is rotated −45° and bored, so it needs the digitizer's rotation and bore
  logic ported over.
- **101.6 mm z-offset root cause** still unknown (§8).

### Code hygiene

- **`prolog.fcl` has duplicated keys.** In the `LaBrPileup` stanza, every key after
  `photonGunMaxZ` is written **twice**. FHiCL takes the last value so it is currently
  harmless, but it should be cleaned up before anyone edits one copy and not the other.
- **`HPGeWaveformsFromStepPointMCs_module_test.cc` is dead.** Its content was copied over the
  real `_module.cc`. It should be deleted to stop it drifting.
- **The digitizers' random engines are never seeded (found 2026-09-29, NOT fixed).**
  `LaBrWaveformsFromStepPointMCsTest` and `HPGeWaveformsFromStepPointMCs` each hold a
  default-constructed `std::default_random_engine _randomGen` that is never given a seed.
  Every job therefore replays the same random sequence for photostatistics / excess noise /
  Fano and ADC noise, whatever `services.SeedService.baseSeed` is set to. Consequences:
  - Changing `baseSeed` changes Geant4 and the generators, but not the digitizer noise sequence.
  - Two jobs with identical deposits give bit-identical waveforms. This is why the old-vs-new
    gun test in `compare.md` matched bin for bin.
  - Grid jobs that differ only in seed share one digitizer noise sequence. Their deposits
    differ, so the results are not simply duplicates, but the smearing is not independent
    between jobs.
  Possible fix, when wanted: seed from a fhicl parameter combined with the run/subrun/event ID
  (reproducible and distinct per event, and avoids the `SeedService` link problem in §12, trap 3).

---

## 12. Working environment and known traps

### Build and run

Edit locally at `~/Desktop/Physics/Mu2e/Simulation/`, copy to the gpvm at
`/exp/mu2e/app/users/sclarke/Simulation/`, build and run there.

```bash
muse setup
muse build -j8

mu2e -c Offline/STMMC/fcl/PoissonPhotonGunToLaBrWaveform.fcl
mu2e -c Offline/STMMC/fcl/PoissonPhotonGunToLaBrWaveform.fcl -n 10000

root -l -b -q 'Offline/STMMC/src/plotLaBrWaveform.C("nts.PoissonPhotonGunToLaBrWaveform.root")'
root -l -b -q 'Offline/STMMC/src/dumpLaBrTruthPulses.C("sim.X.art","truth.csv")'
```

With `verbosityLevel: 1` the digitizer prints all its parameters at startup, including the
derived `peakFraction`, `chargeToADC` and the calibration mode. **Check that printout** —
it is the fastest way to confirm the job is running the configuration you think it is.

Changes to `.fcl` files and `.C` macros need no rebuild. Changes to `.cc` modules do.

### Traps that have cost time

1. **Art plugin file naming.** A module source file **must** be named `ModuleName_module.cc`.
   CMake `cet_build_plugin` only builds those; a `_module_test.cc` is silently **not built**,
   and you keep running the old code. Put the "Test" suffix in the class name, not the
   filename. This one has bitten twice — it is how the HPGe digitizer ran for a long time with
   a crystal centre ~87 mm off in z, rejecting ~99% of steps.

2. **A failed link leaves the old `.so` in place.** If the build errors, the plugin is *not*
   updated, and the job at runtime still shows the **old** module's configuration — for
   example complaining about parameters the new version added. Always confirm the build
   actually **linked**, not merely that the source file is new.

3. **`SeedService` does not link for some STMMC modules.** Adding `Offline::SeedService` to
   CMakeLists does **not** fix it (`undefined reference to mu2e::SeedService::getSeed()`).
   Both digitizers work around it with `std::default_random_engine` plus
   `std::normal_distribution` / `std::poisson_distribution`. `PoissonPhotonGun` and
   `La138Source` *do* use SeedService successfully via `createEngine(...)`. If a new module
   fails to link, do not fight it — use stdlib random, as the digitizers do. But **seed the
   engine yourself**: the digitizers currently do not, so their sequence is the same in every
   job (§11, Code hygiene).

4. **FHiCL is ASCII-only.** Em dashes, smart quotes and similar cause parse errors, including
   inside comments.

5. **Start from the library list of a working module.** When creating a new module by
   adapting an existing one, copy its `cet_build_plugin` `LIBRARIES` block first and only add
   dependencies the code actually needs.

6. **Product instance names.** The digitizers produce *instanced* products — `DigiHPGe:HPGe`
   and `DigiLaBr:LaBr` — not empty instances. Reco prolog tags must use the full instanced
   form. This broke once when the improved HPGe module was promoted.

7. **There is a stale root-level `prolog.fcl` copy** in the working directory that is missing
   the entire `LaBrDigitization` block. **Do not overwrite the in-tree prolog with it** — it
   would break every LaBr job.

8. **A `.bin`/truth mismatch is almost always the per-event timing drift** (§2). 628 samples
   span 1695.6 ns; the timeline advances 1695.0 ns.
