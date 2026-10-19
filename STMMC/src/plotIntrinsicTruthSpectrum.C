// Plot the Geant4 truth energy spectrum of the LaBr crystal split by source: external
// calibration source, La-138 and the Ac-227 chain.
//
// What one histogram entry is:
//   - DEFAULT, from TruthLaBr/decays (LaBrTruthSpectrum analyzer): one DECAY. Its energy is
//     the summed Geant4 deposit of all that decay's products in the crystal (the 1436 keV
//     gamma, X-ray and Auger electron of a La-138 EC decay together) -- what a detector with
//     perfect resolution records per pulse. "total" is all decays together. Pure Geant4
//     truth: no digitizer, no pileup tagging.
//   - FALLBACK, from TruthLaBr/ttree (files made before the per-decay tree existed): one
//     MICROPULSE (1695 ns). Its energy per source is the sum over every decay of that source
//     in the window, so at a boosted activityScale several decays can sum into one entry.
// Written by CalibrationSourcePlusIntrinsicToLaBrWaveform.fcl (or
// CalibrationSourceToLaBrWaveform.fcl, where the Ac-227 histogram is simply empty).
//
// Energies are light-equivalent: alpha deposits are already scaled by alphaQuenchingFactor,
// so the Ac-227 alphas appear where they would in a gamma-calibrated spectrum.
//
// Two pages, log y, 8 keV bins:
//   1) 0-3000 keV, all components stacked as lines, with the key features marked
//   2) 1300-2600 keV zoom: La-138 sum peak, alpha band and the 1808.7 keV signal line
// Counts are raw simulated counts; for physical rates (divided by livetime and
// activityScale) use analyseCalibrationRates.C.
//
// Usage: root -l -b -q 'plotIntrinsicTruthSpectrum.C("nts.CalibrationSourcePlusIntrinsic.root")'
// Original author: Sena Clarke

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"

#include <iostream>

namespace {
  void markLine(double e, double yLow, double yHigh, const char* label) {
    TLine* l = new TLine(e, yLow, e, yHigh);
    l->SetLineStyle(2);
    l->SetLineColor(kGray + 2);
    l->Draw();
    TLatex* t = new TLatex(e, yHigh, label);
    t->SetTextSize(0.025);
    t->SetTextAngle(90);
    t->SetTextAlign(31);
    t->SetTextColor(kGray + 2);
    t->Draw();
  }
}

void plotIntrinsicTruthSpectrum(const char* filename = "nts.CalibrationSourcePlusIntrinsic.root") {
  gStyle->SetOptStat(0);

  TFile* f = TFile::Open(filename);
  if (!f || f->IsZombie()) { std::cerr << "Cannot open " << filename << std::endl; return; }
  TTree* tp = (TTree*)f->Get("TruthLaBr/decays");
  TTree* t  = (TTree*)f->Get("TruthLaBr/ttree");
  const bool perDecay = (tp != nullptr);
  if (!tp && !t) { std::cerr << "No TruthLaBr/decays or TruthLaBr/ttree in " << filename << std::endl; return; }
  std::cout << (perDecay ? "Using TruthLaBr/decays: Geant4 truth, one entry per decay"
                         : "Using TruthLaBr/ttree: Geant4 truth, one entry per micropulse (old file)") << std::endl;

  const int nBins = 375; const double eMax = 3000.0;
  TH1D* hTot = new TH1D("hTot", perDecay ? "LaBr truth spectrum by source (per decay);Deposited energy [keV];Decays / 8 keV"
                                         : "LaBr truth spectrum by source (per micropulse);Deposited energy [keV];Micropulses / 8 keV",
                        nBins, 0, eMax);
  TH1D* hSrc = new TH1D("hSrc", "", nBins, 0, eMax);
  TH1D* hLa  = new TH1D("hLa",  "", nBins, 0, eMax);
  TH1D* hAc  = new TH1D("hAc",  "", nBins, 0, eMax);

  long long nSrc = 0, nLa = 0, nAc = 0;
  if (perDecay) {
    double energy = 0; int source = -1;
    tp->SetBranchAddress("energy", &energy);
    tp->SetBranchAddress("source", &source);
    for (Long64_t i = 0; i < tp->GetEntries(); ++i) {
      tp->GetEntry(i);
      if (energy <= 0) continue;
      hTot->Fill(energy);
      if (source == 0)      { hSrc->Fill(energy); ++nSrc; }
      else if (source == 1) { hLa->Fill(energy);  ++nLa; }
      else if (source == 2) { hAc->Fill(energy);  ++nAc; }
    }
  } else {
    double edepTotal = 0, edepGun = 0, edepLa138 = 0, edepAc227 = 0;
    t->SetBranchAddress("edepTotal", &edepTotal);
    t->SetBranchAddress("edepGun", &edepGun);
    t->SetBranchAddress("edepLa138", &edepLa138);
    const bool haveAc = t->GetBranch("edepAc227") != nullptr;
    if (haveAc) t->SetBranchAddress("edepAc227", &edepAc227);
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (edepTotal > 0) hTot->Fill(edepTotal);
      if (edepGun   > 0) { hSrc->Fill(edepGun);   ++nSrc; }
      if (edepLa138 > 0) { hLa->Fill(edepLa138);  ++nLa; }
      if (haveAc && edepAc227 > 0) { hAc->Fill(edepAc227); ++nAc; }
    }
  }
  std::cout << "entries with deposits -- source: " << nSrc << ", La-138: " << nLa
            << ", Ac-227: " << nAc << std::endl;

  hTot->SetLineColor(kBlack);      hTot->SetLineWidth(2);
  hSrc->SetLineColor(kAzure + 2);  hSrc->SetLineWidth(2);
  hLa->SetLineColor(kOrange + 7);  hLa->SetLineWidth(2);
  hAc->SetLineColor(kGreen + 2);   hAc->SetLineWidth(2);

  auto drawAll = [&](TCanvas* c, double xLo, double xHi) {
    c->SetLogy();
    hTot->GetXaxis()->SetRangeUser(xLo, xHi);
    hTot->SetMinimum(0.5);
    hTot->SetMaximum(hTot->GetMaximum() * 20);
    hTot->Draw("hist");
    hSrc->Draw("hist same");
    hLa->Draw("hist same");
    if (nAc) hAc->Draw("hist same");
    TLegend* leg = new TLegend(0.62, 0.72, 0.89, 0.89);
    leg->AddEntry(hTot, "total", "l");
    leg->AddEntry(hSrc, Form("external source (%lld)", nSrc), "l");
    leg->AddEntry(hLa,  Form("La-138 (%lld)", nLa), "l");
    if (nAc) leg->AddEntry(hAc, Form("Ac-227 chain (%lld)", nAc), "l");
    leg->SetBorderSize(0);
    leg->Draw();
  };

  // Page 1: full range
  TCanvas* c1 = new TCanvas("c1", "full", 1100, 700);
  drawAll(c1, 0, eMax);
  const double yTop = hTot->GetMaximum() * 0.3;
  markLine(37.4, 0.5, yTop, "37 Ba X");
  markLine(788.7, 0.5, yTop, "789");
  markLine(1052.4, 0.5, yTop, "1052 beta end");
  markLine(1473.2, 0.5, yTop, "1473 EC-K");
  markLine(1808.7, 0.5, yTop, "1808.7 signal");
  c1->SaveAs("LaBr_IntrinsicTruth_Full.pdf");
  c1->SaveAs("LaBr_IntrinsicTruth_Full.png");

  // Page 2: signal region zoom
  TCanvas* c2 = new TCanvas("c2", "zoom", 1100, 700);
  drawAll(c2, 1300, 2600);
  const double yTop2 = hTot->GetMaximum() * 0.3;
  markLine(1441.8, 0.5, yTop2, "1442 EC-L");
  markLine(1473.2, 0.5, yTop2, "1473 EC-K");
  markLine(1808.7, 0.5, yTop2, "1808.7 signal");
  c2->SaveAs("LaBr_IntrinsicTruth_SignalRegion.pdf");
  c2->SaveAs("LaBr_IntrinsicTruth_SignalRegion.png");

  std::cout << "wrote LaBr_IntrinsicTruth_Full.{pdf,png} and LaBr_IntrinsicTruth_SignalRegion.{pdf,png}" << std::endl;
}
