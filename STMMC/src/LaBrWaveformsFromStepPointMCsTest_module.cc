// Simulates the electronics response of the LaBr3:Ce+Sr detector.
// Models scintillation light yield, PMT photoelectron conversion with
// dynode-chain excess noise, the single-photoelectron pulse shape
// (1-exp(-t/tau1))^p * exp(-t/tau2) fitted to real LaBr3 data, direct PMT
// readout, noise, and ADC digitization. Generates one STMWaveformDigi per
// micropulse.
//
// Detector parameters follow the Mu2e ELBE test-beam qualification of the
// LaBr3:Ce+Sr / Hamamatsu R10233-100 detector (Huang et al., NIM A draft):
//   - LaBr3:Ce+Sr scintillation decay ~25 ns (slower than standard 16 ns
//     LaBr3:Ce; the Sr co-doping trades speed for energy resolution)
//   - super-bialkali photocathode quantum efficiency ~35%
//   - PMT excess noise factor calibrated so photostatistics reproduce the
//     measured 2.2% FWHM at 662 keV
// Original author: Sena Clarke
// Adapted from HPGe digitization code by Pawel Plesniak

// stdlib includes
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <string>

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

// exception handling
#include "cetlib_except/exception.h"

// fhicl includes
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline includes
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"
#include "Offline/STMMC/inc/LaBrTruthHelpers.hh"

// ROOT includes
#include "art_root_io/TFileService.h"
#include "TTree.h"

namespace mu2e {
  class LaBrWaveformsFromStepPointMCsTest : public art::EDProducer {
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config {
      fhicl::OptionalAtom<art::InputTag> StepPointMCsTagEle{ Name("StepPointMCsTagEle"), Comment("InputTag for StepPointMCs derived from EleBeamCat")};
      fhicl::OptionalAtom<art::InputTag> StepPointMCsTagMu{ Name("StepPointMCsTagMu"), Comment("InputTag for StepPointMCs derived from MuBeamCat")};
      fhicl::OptionalAtom<art::InputTag> StepPointMCsTag1809{ Name("StepPointMCsTag1809"), Comment("InputTag for StepPointMCs derived from TargetStopsCat")};
      fhicl::OptionalAtom<art::InputTag> StepPointMCsTag{ Name("StepPointMCsTag"), Comment("InputTag for StepPointMCs from an arbitrary source (e.g. photon gun)")};

      fhicl::Atom<double> fADC{ Name("fADC"), Comment("ADC operating frequency [MHz]")};
      fhicl::Atom<double> ADCToEnergy{ Name("EnergyPerADCBin"), Comment("ADC energy calibration [keV/bin], applied to the pulse integral")};
      fhicl::OptionalAtom<double> peakKeVPerADC{ Name("peakKeVPerADC"), Comment("If set, calibrate the gain so the pulse PEAK gives energy = peakKeVPerADC * ADC_peak [keV/ADC]. The digitizer computes the peak fraction from its own pulse shape and sampling, so this stays valid if fADC or the pulse shape change. Overrides EnergyPerADCBin.")};
      fhicl::Atom<double> noiseSD{ Name("NoiseSD"), Comment("Standard deviation of ADC noise [mV]. Set this to 0.0 for the ideal case.")};
      fhicl::Atom<double> risingEdgeDecayConstant{ Name("risingEdgeDecayConstant"), Comment("Preamplifier decay time [us]")};
      fhicl::OptionalAtom<double> pulseTau1{ Name("pulseTau1"), Comment("Pulse rise time constant [ns], fitted to real LaBr3 data. Default 22.3543.")};
      fhicl::OptionalAtom<double> pulseTau2{ Name("pulseTau2"), Comment("Pulse decay time constant [ns], fitted to real LaBr3 data. Default 28.0007.")};
      fhicl::OptionalAtom<double> pulseRiseSharpening{ Name("pulseRiseSharpening"), Comment("Pulse rise sharpening exponent p, fitted to real LaBr3 data. Default 6.6596.")};
      fhicl::OptionalAtom<double> lightYield{ Name("lightYield"), Comment("Scintillation photons per MeV. Default 63000.")};
      fhicl::OptionalAtom<double> quantumEfficiency{ Name("quantumEfficiency"), Comment("PMT super-bialkali quantum efficiency (0-1). Default 0.35.")};
      fhicl::OptionalAtom<double> excessNoiseFactor{ Name("excessNoiseFactor"), Comment("PMT excess noise factor F (>=1) from dynode multiplication. Default 1.27 reproduces 2.2% FWHM at 662 keV.")};
      fhicl::OptionalAtom<double> alphaQuenchingFactor{ Name("alphaQuenchingFactor"), Comment("Light yield of alphas relative to electrons at 6 MeV. Default 1 (no quenching); 0.459 for the Mu2e crystal (bench, prolog.fcl), ~0.3 in the literature. Needed with LaBrIntrinsicSource.")};
      fhicl::OptionalAtom<double> alphaQuenchingSlope{ Name("alphaQuenchingSlope"), Comment("Change of the alpha quenching factor per MeV of alpha energy [1/MeV]. Default 0.")};
      fhicl::OptionalAtom<int> microspillBufferLengthCount{ Name("microspillBufferLengthCount"), Comment("Number of microspills to buffer ahead for")};
      fhicl::OptionalAtom<bool> makeTTree{ Name("makeTTree"), Comment("Controls whether to make the TTree")};
      fhicl::OptionalAtom<double> timeOffset{ Name("timeOffset"), Comment("For debugging, adds a time offset [ns]")};
      fhicl::OptionalAtom<bool> foldTime{ Name("foldTime"), Comment("If true, fold step times into one micropulse with fmod(t, micropulseTime) before the buffer cut, as the tracker, calorimeter and CRV digitizers do. Keeps late hits (e.g. slow neutrons) that would otherwise be rejected. Default false.")};
      fhicl::OptionalAtom<uint> resetEventNumber{ Name("resetEventNumber"), Comment("Event ID at which to reset inter-event charge to zero")};
      fhicl::OptionalAtom<uint> verbosityLevel{ Name("verbosityLevel"), Comment("Controls verbosity")};
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit LaBrWaveformsFromStepPointMCsTest(const Parameters& conf);
  private:
    void produce(art::Event& event) override;
    void beginJob();
    void endJob() override;
    void depositCharge(const StepPointMC& step);
    void decayCharge();
    void addNoise();
    void digitize();
    double pulseShape(double dt) const; // single-pe pulse shape, unnormalised
    void buildPulseCDFTable();
    double pulseCDF(double T) const;  // fraction of pe signal arrived by time T

    // fhicl variables
    std::vector<art::InputTag> StepPointMCsTags;
    double fADC = 0;
    double ADCToEnergy = 0;
    double noiseSD = 0;
    double risingEdgeDecayConstant = 0;
    double pulseTau1 = 22.3543;           // pulse rise time constant [ns], fitted to data
    double pulseTau2 = 28.0007;           // pulse decay time constant [ns], fitted to data
    double pulseP = 6.6596;               // pulse rise sharpening exponent, fitted to data
    double peakKeVPerADC = 0.0;           // if >0, peak-based gain calibration [keV/ADC]
    bool usePeakCalib = false;            // derived: peakKeVPerADC > 0
    double peakFraction = 0.0;            // share of a deposit's signal in its tallest bin
    double lightYield = 63000.0;          // scintillation photons per MeV
    double quantumEfficiency = 0.35;      // PMT super-bialkali quantum efficiency
    double excessNoiseFactor = 1.27;      // PMT dynode-chain excess noise factor F
    double alphaQF = 1.0;                 // alpha/electron light ratio at 6 MeV
    double alphaQFSlope = 0.0;            // its change per MeV of alpha energy
    bool makeTTree = false;
    double timeOffset = 0.0;
    bool foldTime = false;                // fold step times into one micropulse
    uint resetEventNumber = 0;
    int verbosityLevel = 0;

    // Derived scintillator constants
    double photoelectronsPerMeV = 0;      // lightYield * quantumEfficiency

    // Experiment constants
    const double feedbackCapacitance = 1e-12;   // [Farads]
    const double micropulseTime = 1695.0;       // [ns]
    const double _e = 1.602176634e-19;          // elementary charge [Coulombs]

    // ADC variables
    double chargeToADC = 0;
    uint nADCs = 0;
    const int16_t ADCMax = static_cast<int16_t>((-1 * std::pow(2, 15)) + 1);
    double ADC = 0;
    uint32_t eventTimeBuffer = 0;

    // LaBr crystal geometry [mm]
    // No rotation (identity), crystal axis along +z.
    // Z corrected from STMGeometryDump (40993.1) by -101.6 mm to match
    // actual G4 step positions, same correction as HPGe (40714.3 -> 40612.7).
    const double crystalCentreX = -3863.3600;
    const double crystalCentreY = 0.0;
    const double crystalCentreZ = 40891.5000;
    CLHEP::Hep3Vector crystalCentrePosition;
    CLHEP::Hep3Vector hitPosition;
    const double crystalL = 76.2;
    const double crystalR = 38.1;
    const double stepPositionTolerance = 0.1;
    const double maxR = crystalR + stepPositionTolerance;
    const double maxZ = crystalL + stepPositionTolerance;

    // Modelling variables
    double hitR = 0;
    double hitZ = 0;
    uint32_t eventTime = 0;
    double tADC = 0;
    double decayExp = 0;                  // exp(-tADC / preampTau), preamplifier IIR decay
    std::vector<double> _pulseCDFTable;   // tabulated normalised CDF of pulseShape
    double pulseCDFStep = 0.05;           // CDF table time step [ns]
    double lastEventEndDecayedCharge = 0;
    int microspillBufferLengthCount = 0;
    const int defaultMicrospillBufferLengthCount = 2;

    // TTree variables
    TTree* ttree;
    double chargeCollected = 0, chargeDecayed = 0;
    uint eventId = 0;
    uint32_t time = 0;
    int16_t ttreeADC = 0;

    // Diagnostic counters
    uint64_t n_events_seen = 0;
    uint64_t n_steps_total = 0;
    uint64_t n_steps_zero_edep = 0;
    uint64_t n_steps_to_deposit = 0;
    uint64_t n_reject_xcut = 0;
    uint64_t n_reject_timecut = 0;
    uint64_t n_folded = 0;
    uint64_t n_reject_bounds = 0;
    uint64_t n_deposit_ok = 0;

    // Data storage vectors
    std::vector<double> _chargeCollected;
    std::vector<double> _chargeDecayed;
    std::vector<double> _chargeCarryOver;
    std::vector<int16_t> _adcs;

    // Random engine (same approach as original HPGe module)
    std::default_random_engine _randomGen;
    std::normal_distribution<double> _noiseDist{0.0, 1.0};
    std::poisson_distribution<int> _poissonDist{1};
  };

  LaBrWaveformsFromStepPointMCsTest::LaBrWaveformsFromStepPointMCsTest(const Parameters& conf)
    : art::EDProducer{conf},
      fADC(conf().fADC()),
      ADCToEnergy(conf().ADCToEnergy()),
      noiseSD(conf().noiseSD()),
      risingEdgeDecayConstant(conf().risingEdgeDecayConstant()) {
        produces<STMWaveformDigiCollection>("LaBr");

        art::InputTag tag;
        if (conf().StepPointMCsTagEle(tag)) { mayConsume<StepPointMCCollection>(tag); StepPointMCsTags.push_back(tag); }
        if (conf().StepPointMCsTagMu(tag))  { mayConsume<StepPointMCCollection>(tag); StepPointMCsTags.push_back(tag); }
        if (conf().StepPointMCsTag1809(tag)){ mayConsume<StepPointMCCollection>(tag); StepPointMCsTags.push_back(tag); }
        if (conf().StepPointMCsTag(tag))    { mayConsume<StepPointMCCollection>(tag); StepPointMCsTags.push_back(tag); }

        if (defaultMicrospillBufferLengthCount < 2)
          throw cet::exception("RANGE", "defaultMicrospillBufferLengthCount has to be more than 1\n");

        crystalCentrePosition.set(crystalCentreX, crystalCentreY, crystalCentreZ);
        tADC = 1e3 / fADC;

        // Assign optional variables
        microspillBufferLengthCount = conf().microspillBufferLengthCount() ? *(conf().microspillBufferLengthCount()) : defaultMicrospillBufferLengthCount;
        verbosityLevel = conf().verbosityLevel() ? *(conf().verbosityLevel()) : 0;
        pulseTau1 = conf().pulseTau1() ? *(conf().pulseTau1()) : 22.3543;
        pulseTau2 = conf().pulseTau2() ? *(conf().pulseTau2()) : 28.0007;
        pulseP = conf().pulseRiseSharpening() ? *(conf().pulseRiseSharpening()) : 6.6596;
        if (pulseTau1 <= 0.0 || pulseTau2 <= 0.0 || pulseP <= 0.0)
          throw cet::exception("RANGE", "pulseTau1, pulseTau2 and pulseRiseSharpening must be > 0\n");
        buildPulseCDFTable();
        lightYield = conf().lightYield() ? *(conf().lightYield()) : 63000.0;
        quantumEfficiency = conf().quantumEfficiency() ? *(conf().quantumEfficiency()) : 0.35;
        excessNoiseFactor = conf().excessNoiseFactor() ? *(conf().excessNoiseFactor()) : 1.27;
        if (excessNoiseFactor < 1.0)
          throw cet::exception("RANGE", "excessNoiseFactor must be >= 1.0\n");

        alphaQF      = conf().alphaQuenchingFactor() ? *(conf().alphaQuenchingFactor()) : 1.0;
        alphaQFSlope = conf().alphaQuenchingSlope()  ? *(conf().alphaQuenchingSlope())  : 0.0;
        if (alphaQF < 0.0)
          throw cet::exception("RANGE", "alphaQuenchingFactor must be >= 0\n");
        photoelectronsPerMeV = lightYield * quantumEfficiency;

        // Buffer sizes
        double _nADCs = (micropulseTime / tADC) + 1;
        nADCs = (int)_nADCs;
        _chargeCollected.assign(nADCs * microspillBufferLengthCount, 0.0);
        _chargeCarryOver.assign(nADCs * (microspillBufferLengthCount - 1), 0.0);
        _chargeDecayed.assign(nADCs, 0.0);
        _adcs.assign(nADCs, 0);

        // Preamplifier IIR decay constant
        decayExp = exp(-tADC / (risingEdgeDecayConstant * 1e3));

        // Convert noise SD from [mV] to charge-carrier units
        noiseSD *= (1e-3 * feedbackCapacitance / _e);

        // Peak fraction: the largest share of a single deposit's signal that
        // lands in one ADC bin. Deposits start on bin boundaries (integer tStart),
        // so this is the same for every pulse and equals max_k of the per-bin
        // pulseCDF difference. Depends on the pulse shape and sampling, so it is
        // computed here rather than assumed.
        peakFraction = 0.0;
        for (uint k = 0; k * tADC < 40.0 * pulseTau2; k++) {
          double frac = pulseCDF((k + 1) * tADC) - pulseCDF(k * tADC);
          if (frac > peakFraction)
            peakFraction = frac;
        }

        // Charge-to-ADC conversion: work in units of photoelectrons.
        // 1 photoelectron = (1e6 / photoelectronsPerMeV) eV, so the integral-based
        // gain is chargeToADC = (1000 / photoelectronsPerMeV) / EnergyPerADCBin,
        // which makes the pulse INTEGRAL read energy at EnergyPerADCBin keV/count.
        peakKeVPerADC = conf().peakKeVPerADC() ? *(conf().peakKeVPerADC()) : 0.0;
        usePeakCalib = (peakKeVPerADC > 0.0);
        if (usePeakCalib) {
          // Calibrate so the pulse PEAK reads energy at peakKeVPerADC keV/count:
          // peak_pe = N_pe * peakFraction and E[keV] = N_pe * 1000 / pePerMeV,
          // so chargeToADC = (1000 / pePerMeV) / (peakKeVPerADC * peakFraction).
          chargeToADC = (1e3 / photoelectronsPerMeV) / (peakKeVPerADC * peakFraction);
        } else {
          double eV_per_pe = 1e6 / photoelectronsPerMeV;
          chargeToADC = eV_per_pe / (ADCToEnergy * 1e3);
        }

        timeOffset = conf().timeOffset() ? *(conf().timeOffset()) : 0.0;
        foldTime = conf().foldTime() ? *(conf().foldTime()) : false;

        makeTTree = conf().makeTTree() ? *(conf().makeTTree()) : false;
        if (makeTTree) {
          art::ServiceHandle<art::TFileService> tfs;
          ttree = tfs->make<TTree>("ttree", "LaBrWaveformsFromStepPointMCs ttree");
          ttree->Branch("chargeCollected", &chargeCollected, "chargeCollected/D");
          ttree->Branch("chargeDecayed", &chargeDecayed, "chargeDecayed/D");
          ttree->Branch("ADC", &ttreeADC, "ADC/S");
          ttree->Branch("eventId", &eventId, "eventId/i");
          ttree->Branch("time", &time, "time/i");
        }

        resetEventNumber = conf().resetEventNumber() ? *(conf().resetEventNumber()) : 0;
      }

  void LaBrWaveformsFromStepPointMCsTest::beginJob() {
    if (verbosityLevel) {
      std::cout << "STM LaBr digitization parameters" << std::endl;
      std::cout << "\tInput parameters" << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "fADC [MHz]"                      << fADC << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "EnergyPerADCBin [keV/bin]"        << ADCToEnergy << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "NoiseSD [mV]"                     << noiseSD / (1e-3 * feedbackCapacitance / _e) << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "risingEdgeDecayConstant [us]"      << risingEdgeDecayConstant << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "pulseTau1 [ns]"                    << pulseTau1 << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "pulseTau2 [ns]"                    << pulseTau2 << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "pulseRiseSharpening"               << pulseP << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "lightYield [photons/MeV]"          << lightYield << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "quantumEfficiency"                 << quantumEfficiency << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "excessNoiseFactor"                 << excessNoiseFactor << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "alphaQuenchingFactor (at 6 MeV)"   << alphaQF << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "alphaQuenchingSlope [1/MeV]"       << alphaQFSlope << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "photoelectronsPerMeV"              << photoelectronsPerMeV << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "foldTime"                          << foldTime << std::endl;
      std::cout << "\tDerived parameters:" << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "tADC [ns]"                        << tADC << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "nADCs"                            << nADCs << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "NoiseSD [charge carriers]"         << noiseSD << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "chargeToADC [bin/pe]"              << chargeToADC << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "peakFraction (signal in tallest bin)" << peakFraction << std::endl;
      if (usePeakCalib) {
        std::cout << std::left << "\t\t" << std::setw(60) << "calibration mode"                 << "PEAK" << std::endl;
        std::cout << std::left << "\t\t" << std::setw(60) << "peakKeVPerADC [keV/ADC]"           << peakKeVPerADC << std::endl;
        std::cout << std::left << "\t\t" << std::setw(60) << "predicted 1809 keV peak [ADC]"     << 1809.0 / peakKeVPerADC << std::endl;
      } else {
        std::cout << std::left << "\t\t" << std::setw(60) << "calibration mode"                 << "INTEGRAL" << std::endl;
        std::cout << std::left << "\t\t" << std::setw(60) << "EnergyPerADCBin [keV/bin]"         << ADCToEnergy << std::endl;
      }
      std::cout << std::left << "\t\t" << std::setw(60) << "pulseCDF table size"               << _pulseCDFTable.size() << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "preampDecayExp (per bin)"          << decayExp << std::endl;
      std::cout << std::endl;
    }
  }

  void LaBrWaveformsFromStepPointMCsTest::produce(art::Event& event) {
    eventId = event.id().event();
    ++n_events_seen;

    for (const auto& tag : StepPointMCsTags) {
      art::Handle<StepPointMCCollection> handle;
      event.getByLabel(tag, handle);
      if (!handle.isValid()) continue;
      for (const StepPointMC& step : *handle) {
        ++n_steps_total;
        if (step.ionizingEdep() != 0) {
          ++n_steps_to_deposit;
          depositCharge(step);
        } else {
          ++n_steps_zero_edep;
        }
      }
    }

    decayCharge();

    if (resetEventNumber != 0 && eventId == resetEventNumber) {
      std::fill(_chargeCarryOver.begin(), _chargeCarryOver.end(), 0);
    }

    addNoise();
    digitize();

    for (auto j : _adcs) {
      if (j > 5000) {
        std::cout << "LaBrDigi WARNING: ADC value " << j << " in event " << eventId << std::endl;
        break;
      }
    }

    eventTimeBuffer = eventId % 5;
    if (eventTimeBuffer == 0 || (eventId % 3) == 0)
      eventTime += nADCs + 1;
    else
      eventTime += nADCs;
    STMWaveformDigi _waveformDigi(eventTime, _adcs);
    std::unique_ptr<STMWaveformDigiCollection> outputDigis(new STMWaveformDigiCollection);
    outputDigis->emplace_back(_waveformDigi);

    if (makeTTree) {
      time = _waveformDigi.trigTimeOffset() - 1;
      for (uint i = 0; i < nADCs; i++) {
        chargeCollected = _chargeCollected[i];
        chargeDecayed = _chargeDecayed[i];
        ttreeADC = _adcs[i];
        time++;
        ttree->Fill();
      }
    }

    _chargeCarryOver.clear();
    _chargeCarryOver.assign(_chargeCollected.begin() + nADCs, _chargeCollected.end());
    _chargeCollected.clear();
    _chargeCollected.assign(_chargeCarryOver.begin(), _chargeCarryOver.end());
    _chargeCollected.insert(_chargeCollected.end(), nADCs, 0);
    std::fill(_chargeDecayed.begin(), _chargeDecayed.end(), 0);
    std::fill(_adcs.begin(), _adcs.end(), 0);

    event.put(std::move(outputDigis), "LaBr");
  }

  void LaBrWaveformsFromStepPointMCsTest::depositCharge(const StepPointMC& step) {
    hitPosition = step.position();

    // LaBr is on the +x side of the STM. Reject HPGe-side hits.
    // HPGe crystal centre is at x ~ -3986, LaBr at x ~ -3863.
    // Midpoint is ~ -3925. Accept hits with x > -3925.
    if (hitPosition.x() < -3925) {
      ++n_reject_xcut;
      return;
    }

    // Time of the hit within the buffer. With foldTime, late hits (e.g. slow neutrons from
    // earlier micropulses) are folded into one micropulse, as in the tracker/calo/CRV digitizers.
    double stepTime = step.time();
    if (foldTime && stepTime >= micropulseTime) {
      stepTime = std::fmod(stepTime, micropulseTime);
      ++n_folded;
    }

    if (stepTime > (microspillBufferLengthCount * micropulseTime)) {
      ++n_reject_timecut;
      return;
    }

    // Transform to crystal-local coordinates.
    // LaBr has identity rotation — just subtract crystal centre.
    hitPosition -= crystalCentrePosition;
    // Shift so front face is at z = 0
    hitPosition.setZ(hitPosition.z() + (crystalL / 2));

    hitR = hitPosition.perp();
    hitZ = hitPosition.z();

    // Bounds check — solid cylinder, no bore hole
    if (hitZ < -stepPositionTolerance || hitZ > maxZ || hitR > maxR) {
      ++n_reject_bounds;
      if (n_reject_bounds <= 20) {
        std::cout << "LaBrDigi bounds-reject #" << n_reject_bounds
                  << " world=" << step.position()
                  << " local(R,Z)=(" << hitR << "," << hitZ << ")"
                  << " maxR=" << maxR << " maxZ=" << maxZ << std::endl;
      }
      return;
    }

    // Scintillator model: energy deposit produces photoelectrons
    // with Poisson statistics, then distributed in time following
    // the scintillation exponential decay profile. Alpha steps are
    // scaled by the alpha/electron light ratio (quenching).
    double N_pe_mean = LaBrTruth::lightEquivalentEdep(step, alphaQF, alphaQFSlope) * photoelectronsPerMeV;
    std::poisson_distribution<int> poisson(N_pe_mean);
    int N_pe = poisson(_randomGen);
    if (N_pe <= 0) {
      ++n_deposit_ok;
      return;
    }

    // PMT excess noise factor: dynode-chain multiplication broadens the
    // measured charge beyond Poisson photostatistics. The output charge has
    // relative variance F/N_pe, so add Gaussian excess variance (F-1)*N_pe on
    // top of the Poisson count. With QE=0.35 and F=1.27 this reproduces the
    // 2.2% FWHM at 662 keV measured for the LaBr3:Ce+Sr / R10233-100 detector.
    double signal_pe = static_cast<double>(N_pe);
    if (excessNoiseFactor > 1.0)
      signal_pe += _noiseDist(_randomGen) * std::sqrt((excessNoiseFactor - 1.0) * N_pe);
    if (signal_pe <= 0.0) {
      ++n_deposit_ok;
      return;
    }

    // Distribute the photoelectron signal across time bins following the
    // single-photoelectron pulse shape fitted to real LaBr3 data (see
    // pulseShape/pulseCDF). Negative sign: signal is a negative-going pulse
    // (same convention as HPGe).
    uint tStart = (stepTime + timeOffset) / tADC;
    uint bufferEnd = nADCs * microspillBufferLengthCount;

    // Fraction of the photoelectron signal remaining to be deposited
    double remaining = signal_pe;
    for (uint t = tStart; t < bufferEnd && remaining > 0.5; t++) {
      double dt = (t - tStart) * tADC;
      double fraction = pulseCDF(dt + tADC) - pulseCDF(dt);
      double deposit = signal_pe * fraction;
      _chargeCollected[t] += -deposit;  // negative-going pulse
      remaining -= deposit;
    }

    ++n_deposit_ok;
  }

  void LaBrWaveformsFromStepPointMCsTest::endJob() {
    std::cout << "===== LaBrWaveformsFromStepPointMCsTest loss breakdown =====" << std::endl;
    std::cout << "  events seen                : " << n_events_seen      << std::endl;
    std::cout << "  steps total                : " << n_steps_total      << std::endl;
    std::cout << "  steps with ionizingEdep==0 : " << n_steps_zero_edep  << std::endl;
    std::cout << "  steps reaching depositCharge: " << n_steps_to_deposit << std::endl;
    std::cout << "    rejected: x < -3925 (HPGe): " << n_reject_xcut     << std::endl;
    std::cout << "    folded into micropulse    : " << n_folded           << std::endl;
    std::cout << "    rejected: step.time > buf : " << n_reject_timecut   << std::endl;
    std::cout << "    rejected: out-of-bounds   : " << n_reject_bounds    << std::endl;
    std::cout << "  deposits successful        : " << n_deposit_ok       << std::endl;
    std::cout << "==============================================================" << std::endl;
  }

  double LaBrWaveformsFromStepPointMCsTest::pulseCDF(double T) const {
    // Fraction of a deposit's photoelectron signal that has arrived by time T
    // (>= 0) after the deposit. The per-bin deposit is pulseCDF(t2)-pulseCDF(t1),
    // and the bins sum to 1 (signal is conserved). Linear interpolation in the
    // tabulated CDF of pulseShape.
    if (T <= 0.0)
      return 0.0;
    double x = T / pulseCDFStep;
    size_t i = static_cast<size_t>(x);
    if (i + 1 >= _pulseCDFTable.size())
      return 1.0;
    double f = x - i;
    return _pulseCDFTable[i] * (1.0 - f) + _pulseCDFTable[i + 1] * f;
  }

  double LaBrWaveformsFromStepPointMCsTest::pulseShape(double dt) const {
    // Single-photoelectron pulse shape fitted to real LaBr3 data:
    //   rise(dt)  = 1 - exp(-dt/tau1)
    //   shape(dt) = rise(dt)^p * exp(-dt/tau2)   for dt > 0 and rise > 0
    //   shape(dt) = 0                            otherwise
    if (dt <= 0.0)
      return 0.0;
    double rise = 1.0 - exp(-dt / pulseTau1);
    if (rise <= 0.0)
      return 0.0;
    return pow(rise, pulseP) * exp(-dt / pulseTau2);
  }

  void LaBrWaveformsFromStepPointMCsTest::buildPulseCDFTable() {
    // The fitted shape has no closed-form integral for non-integer p, so
    // tabulate its normalised cumulative integral (trapezoidal rule) on a
    // pulseCDFStep grid. The tail decays as exp(-t/tau2): 40*tau2 leaves a
    // truncated fraction of order exp(-40), negligible.
    const double tMax = 40.0 * pulseTau2;
    const size_t n = static_cast<size_t>(tMax / pulseCDFStep) + 1;
    _pulseCDFTable.assign(n, 0.0);
    double prev = pulseShape(0.0);
    for (size_t i = 1; i < n; i++) {
      double curr = pulseShape(i * pulseCDFStep);
      _pulseCDFTable[i] = _pulseCDFTable[i - 1] + 0.5 * (prev + curr) * pulseCDFStep;
      prev = curr;
    }
    double total = _pulseCDFTable.back();
    if (total <= 0.0)
      throw cet::exception("RANGE", "pulse shape has non-positive integral\n");
    for (auto& v : _pulseCDFTable)
      v /= total;
  }

  void LaBrWaveformsFromStepPointMCsTest::decayCharge() {
    // Direct PMT model: voltage tracks photoelectron arrival rate,
    // no charge-sensitive preamplifier integration.
    for (uint t = 0; t < nADCs; t++)
      _chargeDecayed[t] = _chargeCollected[t];
  }

  void LaBrWaveformsFromStepPointMCsTest::addNoise() {
    if (noiseSD < std::numeric_limits<double>::epsilon())
      return;
    std::normal_distribution<double> noiseDist(0.0, noiseSD);
    for (size_t i = 0; i < nADCs; i++)
      _chargeDecayed[i] += noiseDist(_randomGen);
  }

  void LaBrWaveformsFromStepPointMCsTest::digitize() {
    for (uint i = 0; i < nADCs; i++) {
      ADC = _chargeDecayed[i] * chargeToADC;
      _adcs[i] = ADC > ADCMax ? static_cast<int16_t>(std::round(ADC)) : ADCMax;
    }
  }
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::LaBrWaveformsFromStepPointMCsTest)
