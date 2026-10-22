// Plot the digitized ADC waveforms of both STM detectors, one pad per event with a pulse.
// Reads the per-ADC-sample TTrees written by the digitizers (makeTTree : true):
//   DigiHPGe/ttree and DigiLaBr/ttree, branches eventId, time [ADC ticks], ADC.
// An event is plotted when its waveform swing (max - min ADC) exceeds minSwing, so noise-only
// events are skipped. Pulses are negative-going. Note: with microspillBufferLengthCount 2, a pulse
// tail can continue into the next event, which then also shows up here.
// Output: STMWaveforms_HPGe.pdf, STMWaveforms_LaBr.pdf (one page per 6 events) + a printed list.
// Usage: root -l -b -q 'plotSTMWaveforms.C("nts.VDResamplerCrystalsHPGeLaBr.root")'
//        root -l -b -q 'plotSTMWaveforms.C("nts.VDResamplerCrystalsHPGeLaBr.root", 50, 12)'
// Original author: Sena Clarke

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TTree.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <vector>

namespace {
  struct Sample { UInt_t time; Short_t adc; };

  void plotDetector(TFile* f, const char* dir, const char* label, double fADC_MHz,
                    double minSwing, int maxEvents) {
    TTree* t = (TTree*)f->Get(Form("%s/ttree", dir));
    if (!t) { std::cout << label << ": no " << dir << "/ttree in file, skipped" << std::endl; return; }

    UInt_t eventId = 0, time = 0;
    Short_t adc = 0;
    t->SetBranchAddress("eventId", &eventId);
    t->SetBranchAddress("time", &time);
    t->SetBranchAddress("ADC", &adc);

    std::map<UInt_t, std::vector<Sample>> events;
    for (Long64_t i = 0; i < t->GetEntries(); i++) {
      t->GetEntry(i);
      events[eventId].push_back({time, adc});
    }

    std::vector<UInt_t> withPulse;
    for (auto& [eid, s] : events) {
      auto mm = std::minmax_element(s.begin(), s.end(),
                                    [](const Sample& a, const Sample& b) { return a.adc < b.adc; });
      double swing = mm.second->adc - mm.first->adc;
      if (swing > minSwing) {
        withPulse.push_back(eid);
        std::cout << label << " event " << eid << ": min ADC " << mm.first->adc
                  << ", max ADC " << mm.second->adc << ", " << s.size() << " samples" << std::endl;
      }
    }
    std::cout << label << ": " << events.size() << " events, " << withPulse.size()
              << " with swing > " << minSwing << " ADC" << std::endl;
    if (withPulse.empty()) return;

    const double tADC = 1e3 / fADC_MHz; // ns per sample
    const int nPlot = std::min<int>(withPulse.size(), maxEvents);
    const TString pdf = Form("STMWaveforms_%s.pdf", label);
    TCanvas* c = new TCanvas(Form("c_%s", label), label, 1200, 800);
    c->Print(pdf + "[");
    for (int p = 0; p < nPlot; p++) {
      if (p % 6 == 0) { c->Clear(); c->Divide(3, 2); }
      auto& s = events[withPulse[p]];
      std::sort(s.begin(), s.end(), [](const Sample& a, const Sample& b) { return a.time < b.time; });
      TGraph* g = new TGraph(s.size());
      for (size_t i = 0; i < s.size(); i++) g->SetPoint(i, (s[i].time - s.front().time) * tADC, s[i].adc);
      g->SetTitle(Form("%s event %u;Time in micropulse [ns];ADC", label, withPulse[p]));
      c->cd(p % 6 + 1);
      g->Draw("AL");
      if (p % 6 == 5 || p == nPlot - 1) c->Print(pdf);
    }
    c->Print(pdf + "]");
    std::cout << label << ": wrote " << pdf << " (" << nPlot << " events)" << std::endl;
  }
}

void plotSTMWaveforms(const char* filename = "nts.VDResamplerCrystalsHPGeLaBr.root",
                      double minSwing = 20, int maxEvents = 60) {
  TFile* f = TFile::Open(filename);
  if (!f || f->IsZombie()) { std::cout << "Cannot open " << filename << std::endl; return; }
  plotDetector(f, "DigiHPGe", "HPGe", 320.0, minSwing, maxEvents);   // STMDAQParameters HPGe fADC
  plotDetector(f, "DigiLaBr", "LaBr", 370.370370370, minSwing, maxEvents); // LaBr fADC
}
