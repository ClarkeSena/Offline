// Simulates the electronic response of the LaBr detector from the geant4 energy deposits
// Original author: Sena Clarke
// Code largely adapted from the HPGe Waveform digitisation code from Pawel Plesniak

// stdlib includes
#include <algorithm>
#include <bits/stdc++.h>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <utility>

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
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
// #include "Offline/DataProducts/inc/STMChannel.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
// #include "Offline/Mu2eUtilities/inc/STMUtils.hh"
// #include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"
// #include "Offline/STMConditions/inc/STMEnergyCalib.hh"

// ROOT includes
#include "art_root_io/TFileService.h"
#include "TTree.h"

namespace mu2e {
    class LaBrWaveformsFromStepPointMCs : public art::EDproducer {
    public:
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        struct Config {
        fhicl::Atom<art::InputTag> StepPointMCsTagEle{ Name("StepPointMCsTagEle"), Comment("InputTag for StepPointMCs derived from EleBeamCat")};
        fhicl::Atom<art::InputTag> StepPointMCsTagMu{ Name("StepPointMCsTagMu"), Comment("InputTag for StepPointMCs derived from MuBeamCat")};
        fhicl::Atom<art::InputTag> StepPointMCsTag1809{ Name("StepPointMCsTag1809"), Comment("InputTag for StepPointMCs derived from TargetStopsCat")};
        fhicl::Atom<double> fADC{ Name("fADC"), Comment("ADC operating frequency [MHz}")};
        fhicl::Atom<double> ADCToEnergy {Name("EnergyPerADCBin"), Comment("ADC energy calibration [keV/bin]")};
        fhicl::Atom<double> noiseSD {Name("NoiseSD"), Comment("Standard deviation of ADC noise [mV]. Set this to 0.0 for the ideal case.")};
        fhicl::Atom<double> risingEdgeDecayConstant{ Name("risingEdgeDecayConstant"), Comment("Rising edge decay time [us]")};
        fhicl::OptionalAtom<int> microspillBufferLengthCount{ Name("microspillBufferLengthCount"), Comment("Number of microspills to buffer ahead for, in number of microspills")};
        fhicl::OptionalAtom<bool> makeTTree{ Name("makeTTree"), Comment("Controls whether to make the TTree with branches chargeCollected, chargeDecayed, ADC, eventId, time")};
        fhicl::OptionalAtom<double> timeOffset{ Name("timeOffset"), Comment("For debugging, adds the named time offset in [ns], used for testing analysis algorithms")};
        fhicl::OptionalAtom<uint> resetEventNumber{ Name("resetEventNumber"), Comment("Simulates the off-spill period by resetting the inter-event last decayed charge to zero")};
        fhicl::OptionalAtom<uint> verbosityLevel{ Name("verbosityLevel"), Comment("Controls verbosity")};
        };
        using Parameters = art::EDProducer::Table<Config>;
        explicit LaBrWaveformsFromStepPointMCs(const Parameters& conf);
    private:
        // fhicl variables
        art::ProductToken<StepPointMCCollection> StepPointMCsTokenEle, StepPointMCsTokenMu, StepPointMCsToken1809;  // Token of StepPointMCs in STMDet
        double fADC = 0;                                                                                            // ADC sampling frequency [MHz]
        double ADCToEnergy = 0;                                                                                     // Calibration of bin width to energy [keV/bin]
        double noiseSD = 0;                                                                                         // Standard deviation of ADC noise [mV]
        double risingEdgeDecayConstant = 0;                                                                         // [us]
        bool makeTTree = false;                                                                                     // Controls whether an analysis TTree is made
        double timeOffset = 0.0;                                                                                    // Used for debugging [ns]
        uint resetEventNumber = 0;                                                                                  // Event ID at which to reset the charge amplitude
        int verbosityLevel = 0;                                                                                     // How much output to generate

        // Define experiment specific constants
        const double feedbackCapacitance = 1e-12;                                                                   // [Farads]
        const double epsilonGe = 2.96;                                                                              // Energy required to generate an eh pair in Ge at 77K [eV]
        const double micropulseTime = 1695.0;                                                                       // [ns]

            // Define physics constants
        const double _e = 1.602176634e-19;                                                                          // Electric charge constant [Coulombs]
        const double electronDriftVelocity = 0.08;                                                                  // Apprixmate charged particle drift velocity [mm/ns]
        const double holeDriftVelocity = 0.06;                                                                      // Apprixmate charged particle drift velocity [mm/ns]

            // ADC variables
        double chargeToADC = 0;                                                                                     // Conversion factor from charge built in capacitor to ADC determined voltage, multiply by this value to get from charge built to ADC voltage.
        uint nADCs = 0;                                                                                             // Number of ADC values in an event
        const int16_t ADCMax = static_cast<int16_t>((-1 * std::pow(2, 15)) + 1);                                    // Maximum ADC value, power is 15 not 16 as using int16_t not uint16_t
        double ADC = 0;                                                                                             // iterator variable
        uint32_t eventTimeBuffer = 0;                                                                               // Multiple of event ids to store

        // Define LaBr variables
            // Define Ge crystal properties [mm]
        const double crystalCentreX = -3863.3600;                                                                     // Crystal centre x position [mm]
        const double crystalCentreY = 0;                                                                            // Crystal centre y position [mm]
        const double crystalCentreZ = 40993.1000;                                                                      // Crystal centre z position [mm]
        CLHEP::Hep3Vector crystalCentrePosition;                                                                    // Crystal centre position vector

        CLHEP::Hep3Vector hitPosition;
        const double crystalL = 76.2;                                                                               // Crystal length [mm]
        const double crystalR = 38.1;                                                                               // Crystal radius [mm]

    }
}