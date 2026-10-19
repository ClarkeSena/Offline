// Simulates the electronics response of the LaBr3:Ce detector.
// Models scintillation light yield, PMT photoelectron conversion,
// exponential scintillation decay profile, preamplifier integration,
// noise, and ADC digitization. Generates one STMWaveformDigi per micropulse.
// Original author: Sena Clarke
// Adapted from HPGe digitization code by Pawel Plesniak

// stdlib includes
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandPoissonQ.h"

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
#include "Offline/SeedService/inc/SeedService.hh"

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
      fhicl::Atom<double> ADCToEnergy{ Name("EnergyPerADCBin"), Comment("ADC energy calibration [keV/bin]")};
      fhicl::Atom<double> noiseSD{ Name("NoiseSD"), Comment("Standard deviation of ADC noise [mV]. Set this to 0.0 for the ideal case.")};
      fhicl::Atom<double> risingEdgeDecayConstant{ Name("risingEdgeDecayConstant"), Comment("Preamplifier decay time [us]")};
      fhicl::OptionalAtom<double> scintDecayTime{ Name("scintDecayTime"), Comment("LaBr3:Ce scintillation decay time [ns]. Default 16 ns.")};
      fhicl::OptionalAtom<double> lightYield{ Name("lightYield"), Comment("Scintillation photons per MeV. Default 63000.")};
      fhicl::OptionalAtom<double> quantumEfficiency{ Name("quantumEfficiency"), Comment("PMT quantum efficiency (0-1). Default 0.25.")};
      fhicl::OptionalAtom<int> microspillBufferLengthCount{ Name("microspillBufferLengthCount"), Comment("Number of microspills to buffer ahead for")};
      fhicl::OptionalAtom<bool> makeTTree{ Name("makeTTree"), Comment("Controls whether to make the TTree")};
      fhicl::OptionalAtom<double> timeOffset{ Name("timeOffset"), Comment("For debugging, adds a time offset [ns]")};
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

    // fhicl variables
    std::vector<art::InputTag> StepPointMCsTags;
    double fADC = 0;
    double ADCToEnergy = 0;
    double noiseSD = 0;
    double risingEdgeDecayConstant = 0;
    double scintDecayTime = 16.0;         // LaBr3:Ce dominant decay time [ns]
    double lightYield = 63000.0;          // scintillation photons per MeV
    double quantumEfficiency = 0.25;      // PMT quantum efficiency
    bool makeTTree = false;
    double timeOffset = 0.0;
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

    // LaBr crystal geometry from STMGeometryDump [mm]
    // No rotation (identity) — crystal axis along +z
    const double crystalCentreX = -3863.3600;
    const double crystalCentreY = 0.0;
    const double crystalCentreZ = 40993.1000;
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
    double scintDecayExp = 0;             // exp(-tADC / scintDecayTime), per-bin scint decay
    double decayExp = 0;                  // exp(-tADC / preampTau), preamplifier IIR decay
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
    uint64_t n_reject_bounds = 0;
    uint64_t n_deposit_ok = 0;

    // Data storage vectors
    std::vector<double> _chargeCollected;
    std::vector<double> _chargeDecayed;
    std::vector<double> _chargeCarryOver;
    std::vector<int16_t> _adcs;

    // Random engine
    art::RandomNumberGenerator::base_engine_t& _engine;
    CLHEP::RandGaussQ _noiseGauss;
    CLHEP::RandPoissonQ _poisson;
  };

  LaBrWaveformsFromStepPointMCsTest::LaBrWaveformsFromStepPointMCsTest(const Parameters& conf)
    : art::EDProducer{conf},
      fADC(conf().fADC()),
      ADCToEnergy(conf().ADCToEnergy()),
      noiseSD(conf().noiseSD()),
      risingEdgeDecayConstant(conf().risingEdgeDecayConstant()),
      _engine(createEngine(art::ServiceHandle<SeedService>()->getSeed())),
      _noiseGauss(_engine, 0.0, 1.0),
      _poisson(_engine, 1.0) {
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
        scintDecayTime = conf().scintDecayTime() ? *(conf().scintDecayTime()) : 16.0;
        lightYield = conf().lightYield() ? *(conf().lightYield()) : 63000.0;
        quantumEfficiency = conf().quantumEfficiency() ? *(conf().quantumEfficiency()) : 0.25;

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

        // Scintillation decay per ADC bin
        scintDecayExp = exp(-tADC / scintDecayTime);

        // Convert noise SD from [mV] to charge-carrier units
        noiseSD *= (1e-3 * feedbackCapacitance / _e);

        // Charge-to-ADC conversion: work in units of photoelectrons
        // 1e3 converts keV to eV
        // For LaBr: 1 photoelectron = (1e6 / photoelectronsPerMeV) eV
        // So chargeToADC = eV_per_pe / (ADCToEnergy * 1e3)
        double eV_per_pe = 1e6 / photoelectronsPerMeV;
        chargeToADC = eV_per_pe / (ADCToEnergy * 1e3);

        timeOffset = conf().timeOffset() ? *(conf().timeOffset()) : 0.0;

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
      std::cout << std::left << "\t\t" << std::setw(60) << "scintDecayTime [ns]"               << scintDecayTime << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "lightYield [photons/MeV]"          << lightYield << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "quantumEfficiency"                 << quantumEfficiency << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "photoelectronsPerMeV"              << photoelectronsPerMeV << std::endl;
      std::cout << "\tDerived parameters:" << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "tADC [ns]"                        << tADC << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "nADCs"                            << nADCs << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "NoiseSD [charge carriers]"         << noiseSD << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "chargeToADC [bin/pe]"              << chargeToADC << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "scintDecayExp (per bin)"           << scintDecayExp << std::endl;
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
      lastEventEndDecayedCharge = 0;
      std::fill(_chargeCarryOver.begin(), _chargeCarryOver.end(), 0);
    } else {
      lastEventEndDecayedCharge = _chargeDecayed.back();
    }

    addNoise();
    digitize();

    for (auto j : _adcs) {
      if (j > 1000)
        throw cet::exception("LogicError", "ADC values too high!");
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

    if (step.time() > (microspillBufferLengthCount * micropulseTime)) {
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
    // the scintillation exponential decay profile.
    double N_pe_mean = step.ionizingEdep() * 1e3 * photoelectronsPerMeV / 1e3;
    // ionizingEdep is in MeV, photoelectronsPerMeV is per MeV
    N_pe_mean = step.ionizingEdep() * photoelectronsPerMeV;
    int N_pe = _poisson.fire(N_pe_mean);
    if (N_pe <= 0) {
      ++n_deposit_ok;
      return;
    }

    // Distribute photoelectrons across time bins following exp(-t/scintDecayTime).
    // The scintillation light arrives at the PMT with this temporal profile.
    // Negative sign: signal is a negative-going pulse (same convention as HPGe).
    uint tStart = (step.time() + timeOffset) / tADC;
    uint bufferEnd = nADCs * microspillBufferLengthCount;

    // Fraction of photoelectrons remaining to be deposited
    double remaining = static_cast<double>(N_pe);
    for (uint t = tStart; t < bufferEnd && remaining > 0.5; t++) {
      double dt = (t - tStart) * tADC;
      double fraction = (1.0 - exp(-(dt + tADC) / scintDecayTime)) -
                        (1.0 - exp(-dt / scintDecayTime));
      double deposit = N_pe * fraction;
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
    std::cout << "    rejected: step.time > buf : " << n_reject_timecut   << std::endl;
    std::cout << "    rejected: out-of-bounds   : " << n_reject_bounds    << std::endl;
    std::cout << "  deposits successful        : " << n_deposit_ok       << std::endl;
    std::cout << "==============================================================" << std::endl;
  }

  void LaBrWaveformsFromStepPointMCsTest::decayCharge() {
    _chargeDecayed[0] = lastEventEndDecayedCharge * decayExp + _chargeCollected[0];
    for (uint t = 1; t < nADCs; t++)
      _chargeDecayed[t] = _chargeDecayed[t-1] * decayExp + _chargeCollected[t];
  }

  void LaBrWaveformsFromStepPointMCsTest::addNoise() {
    if (noiseSD < std::numeric_limits<double>::epsilon())
      return;
    for (size_t i = 0; i < nADCs; i++)
      _chargeDecayed[i] += noiseSD * _noiseGauss.fire();
  }

  void LaBrWaveformsFromStepPointMCsTest::digitize() {
    for (uint i = 0; i < nADCs; i++) {
      ADC = _chargeDecayed[i] * chargeToADC;
      _adcs[i] = ADC > ADCMax ? static_cast<int16_t>(std::round(ADC)) : ADCMax;
    }
  }
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::LaBrWaveformsFromStepPointMCsTest)
