// Plot LaBr digitizer output: full waveform, individual pulses, energy spectrum
// Usage: root -l -b -q 'plotLaBrWaveform.C("nts.PhotonGunToLaBrWaveform.root")'
// Original author: Sena Clarke

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"

#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <iostream>

void plotLaBrWaveform(const char* filename = "nts.PhotonGunToLaBrWaveform.root") {
  gStyle->SetOptStat(0);

  TFile* f = TFile::Open(filename);
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open " << filename << std::endl;
    return;
  }

  TTree* t = (TTree*)f->Get("DigiLaBr/ttree");
  if (!t) {
    std::cerr << "Cannot find DigiLaBr/ttree" << std::endl;
    return;
  }

  int16_t adc;
  UInt_t eventId, timeBin;
  double chargeCollected;
  t->SetBranchAddress("ADC", &adc);
  t->SetBranchAddress("eventId", &eventId);
  t->SetBranchAddress("time", &timeBin);
  t->SetBranchAddress("chargeCollected", &chargeCollected);

  Long64_t nEntries = t->GetEntries();

  // First pass: collect all data and find event boundaries
  struct BinData { UInt_t eventId; UInt_t timeBin; int16_t adc; double charge; };
  std::vector<BinData> allData;
  allData.reserve(nEntries);

  std::map<UInt_t, double> pePerEvent; // eventId -> total photoelectrons deposited
  std::set<UInt_t> eventIds;

  // chargeCollected is in photoelectron units (negative = signal)
  // Convert pe -> keV: E[keV] = N_pe / photoelectronsPerMeV * 1000
  double photoelectronsPerMeV = 63000.0 * 0.35; // lightYield * quantumEfficiency (R10233-100 super-bialkali QE)
  double peToKeV = 1000.0 / photoelectronsPerMeV;

  for (Long64_t i = 0; i < nEntries; i++) {
    t->GetEntry(i);
    allData.push_back({eventId, timeBin, adc, chargeCollected});
    pePerEvent[eventId] += fabs(chargeCollected); // accumulate |pe|
    eventIds.insert(eventId);
  }

  std::cout << "Read " << nEntries << " ADC bins from " << eventIds.size() << " events" << std::endl;

  // =========================================================================
  // Plot 1: Full waveform (all events concatenated)
  // =========================================================================
  TCanvas* c1 = new TCanvas("c1", "Full LaBr Waveform", 1200, 500);
  c1->SetMargin(0.11, 0.03, 0.12, 0.08); // left margin wide enough for the y-axis title

  TGraph* gFull = new TGraph();
  // Find the first signal sample and plot 1000 samples around it
  size_t firstSig = 0;
  for (size_t i = 0; i < allData.size(); i++) {
    if (allData[i].adc < -5) { firstSig = i; break; }
  }
  size_t start = firstSig > 100 ? firstSig - 100 : 0; // 100 samples of baseline before
  size_t nPlot = std::min(allData.size() - start, (size_t)1000);
  for (size_t i = 0; i < nPlot; i++) {
    gFull->SetPoint(i, start + i, allData[start + i].adc);
  }
  gFull->SetTitle("LaBr3:Ce Digitized Waveform (first pulse);ADC sample;ADC value [bins]");
  gFull->SetLineColor(kBlue);
  gFull->SetLineWidth(1);
  gFull->Draw("AL");
  gFull->GetXaxis()->SetTitleSize(0.045);
  gFull->GetYaxis()->SetTitleSize(0.045);
  gFull->GetYaxis()->SetTitleOffset(0.95);

  c1->SaveAs("LaBr_FullWaveform.pdf");
  std::cout << "Saved LaBr_FullWaveform.pdf" << std::endl;

  // =========================================================================
  // Plot 2: Five individual pulses (first 5 events with nonzero signal)
  // =========================================================================
  TCanvas* c2 = new TCanvas("c2", "Individual LaBr Pulses", 1400, 900);
  c2->Divide(3, 2, 0.01, 0.02);

  // Group data by eventId
  std::map<UInt_t, std::vector<std::pair<UInt_t, int16_t>>> eventData;
  for (const auto& d : allData) {
    eventData[d.eventId].push_back({d.timeBin, d.adc});
  }

  // Find events with actual signal (min ADC < -5)
  std::vector<UInt_t> signalEvents;
  for (const auto& [eid, bins] : eventData) {
    int16_t minADC = 0;
    for (const auto& [t, a] : bins) {
      if (a < minADC) minADC = a;
    }
    if (minADC < -5) signalEvents.push_back(eid);
    if (signalEvents.size() >= 5) break;
  }

  TGraph* gPulse[5];
  for (int p = 0; p < 5 && p < (int)signalEvents.size(); p++) {
    c2->cd(p + 1);
    gPad->SetMargin(0.17, 0.05, 0.12, 0.08); // left margin wide enough for the y-axis title

    const auto& bins = eventData[signalEvents[p]];
    double adcToNs = 1e3 / 370.370370370; // time bin to ns

    // Find the pulse: first/last bin above threshold (same threshold as event selection)
    double tFirst = -1, tLast = -1;
    for (const auto& [tb, a] : bins) {
      if (a < -5) {
        double tNs = tb * adcToNs;
        if (tFirst < 0) tFirst = tNs;
        tLast = tNs;
      }
    }

    // Padding before/after the pulse, in ns — adjust these to widen the view
    double padBefore = 20.0;
    double padAfter = 100.0;

    // Fill the graph with only the samples in the zoom window, with time relative
    // to pulse onset, so Draw("AL") auto-ranges to exactly the pulse
    gPulse[p] = new TGraph();
    int np = 0;
    for (const auto& [tb, a] : bins) {
      double tNs = tb * adcToNs;
      if (tFirst >= 0 && (tNs < tFirst - padBefore || tNs > tLast + padAfter)) continue;
      gPulse[p]->SetPoint(np++, tNs - tFirst, a);
    }
    gPulse[p]->SetTitle(
        Form("Event %u;Time since pulse onset [ns];ADC value [bins]", signalEvents[p]));
    gPulse[p]->SetLineColor(kBlue);
    gPulse[p]->SetLineWidth(1);
    gPulse[p]->Draw("AL");
    gPulse[p]->GetXaxis()->SetTitleSize(0.05);
    gPulse[p]->GetYaxis()->SetTitleSize(0.05);
    gPulse[p]->GetYaxis()->SetTitleOffset(1.5);
  }

  // Panel 6: legend/info
  c2->cd(6);
  TLatex* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.08);
  tex->DrawLatex(0.1, 0.8, "LaBr3:Ce Digitizer");
  tex->SetTextSize(0.06);
  tex->DrawLatex(0.1, 0.65, "1.809 MeV photon gun");
  tex->DrawLatex(0.1, 0.50, Form("f_{ADC} = 370.37 MHz"));
  tex->DrawLatex(0.1, 0.35, Form("%lu events total", eventIds.size()));

  c2->SaveAs("LaBr_IndividualPulses.pdf");
  std::cout << "Saved LaBr_IndividualPulses.pdf" << std::endl;

  // =========================================================================
  // Plot 3: Energy spectrum from truth (Geant4 deposited energy per event)
  // =========================================================================
  TCanvas* c3 = new TCanvas("c3", "LaBr Energy Spectrum", 800, 600);
  c3->SetMargin(0.12, 0.05, 0.12, 0.08);

  // Convert pe sums to keV
  std::map<UInt_t, double> energyPerEvent;
  for (const auto& [eid, pe] : pePerEvent) {
    energyPerEvent[eid] = pe * peToKeV;
  }

  // Fixed binning matched to plotLa138TruthSpectrum.C so the digitized spectrum can be
  // compared bin-for-bin with the Geant4 truth spectrum.
  const int    nBins   = 250;
  const double histMin = 0.0;
  const double histMax = 2000.0; // keV
  TH1D* hEnergy = new TH1D("hEnergy",
    "LaBr3:Ce Energy Spectrum (digitized);Energy deposited [keV];Counts / bin",
    nBins, histMin, histMax);

  for (const auto& [eid, e] : energyPerEvent) {
    if (e > 0) hEnergy->Fill(e);
  }

  hEnergy->SetLineColor(kBlue);
  hEnergy->SetLineWidth(2);
  hEnergy->Draw();
  hEnergy->GetXaxis()->SetTitleSize(0.045);
  hEnergy->GetYaxis()->SetTitleSize(0.045);

  // Label photopeak
  int peakBin = hEnergy->GetMaximumBin();
  double peakE = hEnergy->GetBinCenter(peakBin);
  TLatex* label = new TLatex();
  label->SetNDC();
  label->SetTextSize(0.035);
  label->DrawLatex(0.55, 0.85, Form("Photopeak: %.0f keV", peakE));
  label->DrawLatex(0.55, 0.78, Form("Events: %d", (int)hEnergy->GetEntries()));

  c3->SaveAs("LaBr_EnergySpectrum.pdf");
  std::cout << "Saved LaBr_EnergySpectrum.pdf" << std::endl;

  std::cout << "\nDone. Output files:" << std::endl;
  std::cout << "  LaBr_FullWaveform.pdf" << std::endl;
  std::cout << "  LaBr_IndividualPulses.pdf" << std::endl;
  std::cout << "  LaBr_EnergySpectrum.pdf" << std::endl;

  f->Close();
}