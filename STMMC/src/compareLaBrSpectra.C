// Overlay two LaBr runs (e.g. an old job fcl vs the modular LaBrPipeline) to check they agree
// within statistics. Both spectra are normalised PER EVENT (micropulse), so runs with
// different event counts compare directly.
//   1. digitized energy per event, from DigiLaBr/ttree chargeCollected (as plotLaBrWaveform.C)
//   2. Geant4 truth energy per decay, split by source, from TruthLaBr/decays (if both have it)
// Prints chi2 / KS probabilities and per-source rates (decays per event).
//
// Usage: root -l -b -q 'compareLaBrSpectra.C("old.root","new.root","Gun")'
//   -> LaBrCompare_<tag>_Digitized.pdf, LaBrCompare_<tag>_Truth_<source>.pdf
// Original author: Sena Clarke

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
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
  const double eMax = 2000.0; // keV, same binning as plotLaBrWaveform.C
  const double peToKeV = 1000.0 / (63000.0 * 0.35);

  // Per-event digitized energy; returns the number of events in the file (0 on failure)
  long digitized(TFile* f, TH1D* h) {
    TTree* t = (TTree*)f->Get("DigiLaBr/ttree");
    if (!t) { std::cerr << "no DigiLaBr/ttree in " << f->GetName() << "\n"; return 0; }
    UInt_t eventId; double charge;
    t->SetBranchStatus("*", 0);
    t->SetBranchStatus("eventId", 1);
    t->SetBranchStatus("chargeCollected", 1);
    t->SetBranchAddress("eventId", &eventId);
    t->SetBranchAddress("chargeCollected", &charge);
    std::map<UInt_t, double> pe;
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      pe[eventId] += std::fabs(charge);
    }
    for (auto const& kv : pe)
      if (kv.second > 0) h->Fill(kv.second * peToKeV);
    return (long)pe.size();
  }

  // Truth energy per decay by source; fills counts[source]
  bool truth(TFile* f, TH1D* h[3], long counts[3]) {
    TTree* t = (TTree*)f->Get("TruthLaBr/decays");
    if (!t) return false;
    double energy; int source;
    t->SetBranchAddress("energy", &energy);
    t->SetBranchAddress("source", &source);
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (source < 0 || source > 2) continue;
      h[source]->Fill(energy);
      ++counts[source];
    }
    return true;
  }

  void style(TH1D* h, int colour, int lineStyle) {
    h->SetLineColor(colour);
    h->SetLineStyle(lineStyle);
    h->SetLineWidth(2);
  }

  void compare(const char* what, TH1D* a, TH1D* b) {
    if (a->GetEntries() < 1 || b->GetEntries() < 1) return;
    std::cout << "  " << what << ": entries old " << a->GetEntries() << ", new " << b->GetEntries()
              << " | chi2 prob " << a->Chi2Test(b, "UU") << " | KS prob " << a->KolmogorovTest(b) << "\n";
  }

  // Overlay (top, log y) + new/old ratio (bottom)
  void overlay(const char* pdf, const char* title, TH1D* a, TH1D* b, double nA, double nB) {
    TCanvas c("c", title, 900, 800);
    TPad top("top", "", 0, 0.3, 1, 1), bot("bot", "", 0, 0, 1, 0.3);
    top.SetBottomMargin(0.02); top.SetLogy(); bot.SetTopMargin(0.02); bot.SetBottomMargin(0.3);
    top.Draw(); bot.Draw();
    TH1D* sa = (TH1D*)a->Clone(); sa->Scale(1.0 / nA);
    TH1D* sb = (TH1D*)b->Clone(); sb->Scale(1.0 / nB);
    top.cd();
    sa->SetTitle(Form("%s;;Counts / event / bin", title));
    sa->SetMinimum(0.5 / std::max(nA, nB));
    sa->Draw("hist");
    sb->Draw("hist same");
    TLegend leg(0.62, 0.75, 0.88, 0.88);
    leg.AddEntry(sa, Form("old (%.0f events)", nA), "l");
    leg.AddEntry(sb, Form("new (%.0f events)", nB), "l");
    leg.Draw();
    bot.cd();
    TH1D* r = (TH1D*)sb->Clone("ratio");
    r->Divide(sa);
    r->SetTitle(";Energy [keV];new / old");
    r->SetMinimum(0.0); r->SetMaximum(2.0);
    r->SetLineColor(kBlack); r->SetLineStyle(1);
    r->GetXaxis()->SetTitleSize(0.12); r->GetXaxis()->SetLabelSize(0.1);
    r->GetYaxis()->SetTitleSize(0.1); r->GetYaxis()->SetLabelSize(0.08); r->GetYaxis()->SetTitleOffset(0.4);
    r->Draw("e");
    TLine one(0, 1, eMax, 1); one.SetLineStyle(2); one.Draw();
    c.SaveAs(pdf);
  }
}

void compareLaBrSpectra(const char* oldFile, const char* newFile, const char* tag = "Compare") {
  gStyle->SetOptStat(0);
  TH1D::SetDefaultSumw2();
  TFile* fo = TFile::Open(oldFile);
  TFile* fn = TFile::Open(newFile);
  if (!fo || fo->IsZombie() || !fn || fn->IsZombie()) { std::cerr << "cannot open inputs\n"; return; }

  TH1D* dO = new TH1D("dO", "", nBins, 0, eMax);
  TH1D* dN = new TH1D("dN", "", nBins, 0, eMax);
  const long nO = digitized(fo, dO), nN = digitized(fn, dN);
  if (!nO || !nN) return;
  style(dO, kBlue, 1); style(dN, kRed, 2);

  std::cout << "\n=== " << tag << ": " << oldFile << " (old) vs " << newFile << " (new)\n";
  std::cout << "  events: old " << nO << ", new " << nN << "\n";
  std::cout << "  events with a deposit: old " << dO->GetEntries() / nO << ", new " << dN->GetEntries() / nN << " per event\n";
  std::cout << "  mean digitized energy of hit events [keV]: old " << dO->GetMean() << ", new " << dN->GetMean() << "\n";
  compare("digitized", dO, dN);
  overlay(Form("LaBrCompare_%s_Digitized.pdf", tag), Form("%s: digitized energy per event", tag), dO, dN, nO, nN);

  const char* names[3] = {"external", "La-138", "Ac-227"};
  TH1D *tO[3], *tN[3];
  long cO[3] = {0, 0, 0}, cN[3] = {0, 0, 0};
  for (int s = 0; s < 3; ++s) {
    tO[s] = new TH1D(Form("tO%d", s), "", nBins, 0, eMax); style(tO[s], kBlue, 1);
    tN[s] = new TH1D(Form("tN%d", s), "", nBins, 0, eMax); style(tN[s], kRed, 2);
  }
  if (truth(fo, tO, cO) && truth(fn, tN, cN)) {
    std::cout << "  truth decays per event (+- stat), by source:\n";
    for (int s = 0; s < 3; ++s) {
      if (!cO[s] && !cN[s]) continue;
      std::cout << "    " << names[s] << ": old " << double(cO[s]) / nO << " +- " << std::sqrt(double(cO[s])) / nO
                << ", new " << double(cN[s]) / nN << " +- " << std::sqrt(double(cN[s])) / nN << "\n";
      compare(Form("truth %s", names[s]), tO[s], tN[s]);
    }
    for (int s = 0; s < 3; ++s) {
      if (!cO[s] || !cN[s]) continue;
      overlay(Form("LaBrCompare_%s_Truth_%s.pdf", tag, s == 0 ? "External" : (s == 1 ? "La138" : "Ac227")),
              Form("%s: truth energy per decay, %s", tag, names[s]), tO[s], tN[s], nO, nN);
    }
  } else
    std::cout << "  (no TruthLaBr/decays in one of the files: truth comparison skipped)\n";
}
