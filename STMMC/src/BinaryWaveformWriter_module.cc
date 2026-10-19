// Writes STMWaveformDigi ADC samples to a flat binary file: little-endian
// int16, no header, waveforms concatenated in event order. This matches the
// raw-stream format consumed by the offline matched-filter pulse finder
// (e.g. labr_sim_20kHz_ADC.bin); read with numpy as np.fromfile(f, dtype="<i2").
// Original author: Sena Clarke

// stdlib includes
#include <fstream>
#include <iostream>

// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// exception handling
#include "cetlib_except/exception.h"

// fhicl includes
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline includes
#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"

namespace mu2e {
  class BinaryWaveformWriter : public art::EDAnalyzer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> stmWaveformDigisTag{Name("stmWaveformDigisTag"), Comment("Tag identifying the STMWaveformDigis")};
        fhicl::Atom<std::string> fileName{Name("fileName"), Comment("Output binary file name")};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit BinaryWaveformWriter(const Parameters& conf);
      void analyze(const art::Event& e);
      void endJob();
    private:
      art::ProductToken<STMWaveformDigiCollection> stmWaveformDigisToken;
      std::string fileName;
      std::ofstream file;
      size_t samplesWritten = 0;
  };

  BinaryWaveformWriter::BinaryWaveformWriter(const Parameters& conf) :
    art::EDAnalyzer(conf),
    stmWaveformDigisToken(consumes<STMWaveformDigiCollection>(conf().stmWaveformDigisTag())),
    fileName(conf().fileName()) {
      file.open(fileName, std::ios::out | std::ios::binary);
      if (!file.is_open())
        throw cet::exception("FILE") << "Cannot open output binary file " << fileName << "\n";
    };

  void BinaryWaveformWriter::analyze(const art::Event& event) {
    auto const& waveformDigis = event.getProduct(stmWaveformDigisToken);
    for (const STMWaveformDigi& waveformDigi : waveformDigis) {
      const std::vector<int16_t>& adcs = waveformDigi.adcs();
      file.write(reinterpret_cast<const char*>(adcs.data()), adcs.size() * sizeof(int16_t));
      samplesWritten += adcs.size();
    }
  };

  void BinaryWaveformWriter::endJob() {
    file.close();
    std::cout << "BinaryWaveformWriter: wrote " << samplesWritten << " int16 samples ("
              << samplesWritten * sizeof(int16_t) << " bytes) to " << fileName << std::endl;
  };
}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::BinaryWaveformWriter)
