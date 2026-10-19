// Compare the merged HPGe + LaBr job (PhotonGunToHPGeAndLaBrWaveforms.fcl) with the two
// standalone jobs (PhotonGunToHPGeWaveform.fcl, LaBrPipeline.fcl), and show what the same beam
// looks like in both detectors. All spectra are normalised PER EVENT (micropulse), so runs
// with different event counts compare directly.
//
// Energy per event = sum of |chargeCollected| over the event's ADC bins, converted to keV
// (same estimator as plotHPGeAndLaBrSpectra.C / plotLaBrWaveform.C, pre-ADC, carries the
// detector smearing but not ADC quantization or electronic noise):
//   HPGe : e-h pairs x 2.96 eV       (DigiHPGe/ttree)
//   LaBr : pe / (63000 x 0.35) MeV   (DigiLaBr/ttree)
//
// Outputs:
//   MergedHPGeLaBr_Spectra.pdf      merged file: HPGe and LaBr for the same beam, overlaid
//   MergedHPGeLaBr_Correlation.pdf  merged file: E_LaBr vs E_HPGe per event
//   MergedHPGeLaBr_HPGe.pdf         standalone HPGe vs merged HPGe, with ratio (if hpgeFile)
//   MergedHPGeLaBr_LaBr.pdf         standalone LaBr vs merged LaBr, with ratio (if labrFile)
// and prints, per detector, the hit fraction, full-energy-peak fraction and chi2 / KS probs.
//
// Usage: root -l -b -q 'compareMergedHPGeLaBr.C("nts.PhotonGunToHPGeAndLaBrWaveforms.root","nts.PhotonGunToHPGeWaveform.root","nts.LaBrPipeline.root")'
//        Pass "" for a standalone file to skip that comparison.
// Original author: Sena Clarke

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"
#include "TTree.h"

#include <cmath>
#include <iostream>
#include <map>

namespace {
  const int nBins = 250;
  const double eMax = 2000.0;  // keV, same binning as plotHPGeAndLaBrSpectra.C
  const double eBeam = 1809.0; // keV, MergedBeam.E

  struct Detector {
    const char* name;
    const char* tree;
    double keVPerCharge;
    double peakHalfWidth; // keV, full-energy-peak window around eBeam
  };
  const Detector HPGe = {"HPGe", "DigiHPGe/ttree", 2.96 / 1e3, 10.0};
  const Detector LaBr = {"LaBr", "DigiLaBr/ttree", 1e3 / (63000.0 * 0.35), 40.0};

  // Per-event energy [keV] keyed by eventId; every event is present (0 if no deposit)
  bool energies(TFile* f, const Detector& d, std::map<UInt_t, double>& e) {
    TTree* t = (TTree*)f->Get(d.tree);
    if (!t) {
      std::cerr << "no " << d.tree << " in " << f->GetName() << "\n";
      return false;
    }
    UInt_t eventId;
    double charge;
    t->SetBranchStatus("*", 0);
    t->SetBranchStatus("eventId", 1);
    t->SetBranchStatus("chargeCollected", 1);
    t->SetBranchAddress("eventId", &eventId);
    t->SetBranchAddress("chargeCollected", &charge);
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      e[eventId] += std::fabs(charge) * d.keVPerCharge;
    }
    return !e.empty();
  }

  TH1D* spectrum(const std::map<UInt_t, double>& e, const char* name, int colour, int lineStyle) {
    TH1D* h = new TH1D(name, "", nBins, 0, eMax);
    for (auto const& kv : e)
      if (kv.second > 0) h->Fill(kv.second);
    h->SetLineColor(colour);
    h->SetLineStyle(lineStyle);
    h->SetLineWidth(2);
    return h;
  }

  void summary(const char* label, const Detector& d, const std::map<UInt_t, double>& e) {
    const double n = e.size();
    long hit = 0, peak = 0;
    for (auto const& kv : e) {
      if (kv.second > 0) ++hit;
      if (std::fabs(kv.second - eBeam) < d.peakHalfWidth) ++peak;
    }
    std::cout << "  " << label << ": " << e.size() << " events, hit " << hit / n << " +- "
              << std::sqrt(double(hit)) / n << ", full-energy peak (+-" << d.peakHalfWidth
              << " keV) " << peak / n << " +- " << std::sqrt(double(peak)) / n << " per event\n";
  }

  // Standalone (a) vs merged (b), per event: overlay (top, log y) + merged/standalone (bottom)
  void overlay(const char* pdf, const Detector& d, TH1D* a, TH1D* b, double nA, double nB) {
    if (a->GetEntries() > 0 && b->GetEntries() > 0)
      std::cout << "  " << d.name << " standalone vs merged: chi2 prob " << a->Chi2Test(b, "UU")
                << " | KS prob " << a->KolmogorovTest(b) << "\n";
    TCanvas c("c", d.name, 900, 800);
    TPad top("top", "", 0, 0.3, 1, 1), bot("bot", "", 0, 0, 1, 0.3);
    top.SetBottomMargin(0.02); top.SetLogy(); bot.SetTopMargin(0.02); bot.SetBottomMargin(0.3);
    top.Draw(); bot.Draw();
    TH1D* sa = (TH1D*)a->Clone(); sa->Scale(1.0 / nA);
    TH1D* sb = (TH1D*)b->Clone(); sb->Scale(1.0 / nB);
    top.cd();
    sa->SetTitle(Form("%s digitized energy per event;;Counts / event / bin", d.name));
    sa->SetMinimum(0.5 / std::max(nA, nB));
    sa->Draw("hist");
    sb->Draw("hist same");
    TLegend leg(0.55, 0.75, 0.88, 0.88);
    leg.AddEntry(sa, Form("standalone (%.0f events)", nA), "l");
    leg.AddEntry(sb, Form("merged (%.0f events)", nB), "l");
    leg.Draw();
    bot.cd();
    TH1D* r = (TH1D*)sb->Clone("ratio");
    r->Divide(sa);
    r->SetTitle(";Energy [keV];merged / standalone");
    r->SetMinimum(0.0); r->SetMaximum(2.0);
    r->SetLineColor(kBlack); r->SetLineStyle(1);
    r->GetXaxis()->SetTitleSize(0.12); r->GetXaxis()->SetLabelSize(0.1);
    r->GetYaxis()->SetTitleSize(0.1); r->GetYaxis()->SetLabelSize(0.08); r->GetYaxis()->SetTitleOffset(0.4);
    r->Draw("e");
    TLine one(0, 1, eMax, 1); one.SetLineStyle(2); one.Draw();
    c.SaveAs(pdf);
  }

  // Standalone file vs the merged detector's spectrum
  void compareStandalone(const char* file, const Detector& d, TH1D* merged, double nMerged) {
    if (!file || !*file) return;
    TFile* f = TFile::Open(file);
    if (!f || f->IsZombie()) { std::cerr << "cannot open " << file << "\n"; return; }
    std::map<UInt_t, double> e;
    if (!energies(f, d, e)) return;
    std::cout << "\n=== " << d.name << ": " << file << " (standalone)\n";
    summary("standalone", d, e);
    TH1D* h = spectrum(e, Form("h%sStandalone", d.name), kBlue, 1);
    overlay(Form("MergedHPGeLaBr_%s.pdf", d.name), d, h, merged, e.size(), nMerged);
  }
}

void compareMergedHPGeLaBr(const char* mergedFile = "nts.PhotonGunToHPGeAndLaBrWaveforms.root",
                           const char* hpgeFile = "", const char* labrFile = "") {
  gStyle->SetOptStat(0);
  TH1D::SetDefaultSumw2();
  TFile* fm = TFile::Open(mergedFile);
  if (!fm || fm->IsZombie()) { std::cerr << "cannot open " << mergedFile << "\n"; return; }

  std::map<UInt_t, double> eH, eL;
  const bool haveH = energies(fm, HPGe, eH), haveL = energies(fm, LaBr, eL);
  if (!haveH && !haveL) return;

  std::cout << "\n=== merged: " << mergedFile << "\n";
  if (haveH) summary("HPGe", HPGe, eH);
  if (haveL) summary("LaBr", LaBr, eL);
  TH1D* hH = haveH ? spectrum(eH, "hHPGeMerged", kRed + 1, 2) : nullptr;
  TH1D* hL = haveL ? spectrum(eL, "hLaBrMerged", kRed + 1, 2) : nullptr;

  // Same beam, both detectors, per event
  if (haveH && haveL) {
    long both = 0;
    TH2D* corr = new TH2D("corr", "Same beam, per event;E_{HPGe} [keV];E_{LaBr} [keV]",
                          nBins / 2, 0, eMax, nBins / 2, 0, eMax);
    for (auto const& kv : eH) {
      const double l = eL.count(kv.first) ? eL.at(kv.first) : 0.0;
      if (kv.second > 0 && l > 0) ++both;
      if (kv.second > 0 || l > 0) corr->Fill(kv.second, l);
    }
    std::cout << "  both detectors hit in the same event: " << double(both) / eH.size() << " per event\n";

    TCanvas c("cSpectra", "", 900, 650);
    c.SetLogy();
    TH1D* sH = (TH1D*)hH->Clone("sH"); sH->Scale(1.0 / eH.size());
    TH1D* sL = (TH1D*)hL->Clone("sL"); sL->Scale(1.0 / eL.size());
    sH->SetLineColor(kRed + 1); sH->SetLineStyle(1);
    sL->SetLineColor(kBlue + 1); sL->SetLineStyle(1);
    sL->SetTitle("Same beam in both detectors;Energy deposited [keV];Counts / event / bin");
    sL->SetMaximum(2.0 * std::max(sH->GetMaximum(), sL->GetMaximum()));
    sL->SetMinimum(0.5 / std::max(eH.size(), eL.size()));
    sL->Draw("hist");
    sH->Draw("hist same");
    TLegend leg(0.15, 0.75, 0.5, 0.88);
    leg.SetBorderSize(0);
    leg.AddEntry(sL, "LaBr3:Ce (scintillator)", "l");
    leg.AddEntry(sH, "HPGe (semiconductor)", "l");
    leg.Draw();
    c.SaveAs("MergedHPGeLaBr_Spectra.pdf");

    TCanvas c2("cCorr", "", 800, 700);
    c2.SetLogz();
    c2.SetRightMargin(0.14);
    corr->Draw("colz");
    c2.SaveAs("MergedHPGeLaBr_Correlation.pdf");
  }

  if (haveH) compareStandalone(hpgeFile, HPGe, hH, eH.size());
  if (haveL) compareStandalone(labrFile, LaBr, hL, eL.size());
}
