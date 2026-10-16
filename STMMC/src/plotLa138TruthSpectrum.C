// Plot the Geant4 truth energy spectrum of the PhotonGunPlusLa138 simulation, split by
// source. The LaBrTruthSpectrum analyzer classifies each event's LaBr-crystal energy into
// photon-gun signal (1.809 MeV gamma) and La-138 intrinsic activity, by tracing every step
// to its primary. This macro makes three energy spectra from that truth tree:
//   1) combined (all energy deposited in the crystal)
//   2) photon-gun-only (the 1809 keV peak + Compton/escape structure)
//   3) La-138-only (1468/1436 keV EC sum, 789 keV, beta continuum, X-rays)
// Usage: root -l -b -q 'plotLa138TruthSpectrum.C("nts.PhotonGunPlusLa138.root")'
// Original author: Sena Clarke

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"

#include <iostream>

void plotLa138TruthSpectrum(const char* filename = "nts.PhotonGunPlusLa138.root") {
  gStyle->SetOptStat(0);

  TFile* f = TFile::Open(filename);
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open " << filename << std::endl;
    return;
  }

  TTree* t = (TTree*)f->Get("TruthLaBr/ttree");
  if (!t) {
    std::cerr << "Cannot find TruthLaBr/ttree (run PhotonGunPlusLa138ToLaBrWaveform.fcl "
                 "with the LaBrTruthSpectrum analyzer enabled)" << std::endl;
    return;
  }

  double edepTotal, edepGun, edepLa138;
  t->SetBranchAddress("edepTotal", &edepTotal);
  t->SetBranchAddress("edepGun", &edepGun);
  t->SetBranchAddress("edepLa138", &edepLa138);

  Long64_t nEntries = t->GetEntries();

  // Fixed binning over a range that covers the 1809 keV signal photopeak with headroom.
  const int    nBins   = 250;
  const double histMin = 0.0;
  const double histMax = 2000.0; // keV

  TH1D* hComb  = new TH1D("hComb",
    "LaBr3:Ce Truth Spectrum (combined);Energy deposited [keV];Counts / bin",
    nBins, histMin, histMax);
  TH1D* hGun   = new TH1D("hGun",
    "LaBr3:Ce Truth Spectrum (photon gun);Energy deposited [keV];Counts / bin",
    nBins, histMin, histMax);
  TH1D* hLa138 = new TH1D("hLa138",
    "LaBr3:Ce Truth Spectrum (La-138);Energy deposited [keV];Counts / bin",
    nBins, histMin, histMax);

  long nGunEvts = 0, nLa138Evts = 0;
  for (Long64_t i = 0; i < nEntries; i++) {
    t->GetEntry(i);
    // Each generator deposits energy independently; fill its own histogram only when it
    // actually deposited something in this event (the gun fires every event, La-138 rarely).
    if (edepTotal > 0) hComb->Fill(edepTotal);
    if (edepGun   > 0) { hGun->Fill(edepGun);     ++nGunEvts;   }
    if (edepLa138 > 0) { hLa138->Fill(edepLa138); ++nLa138Evts; }
  }

  std::cout << "Read " << nEntries << " events: "
            << nGunEvts << " with photon-gun energy, "
            << nLa138Evts << " with La-138 energy" << std::endl;

  // ---- Plot 1: combined -----------------------------------------------------
  TCanvas* c1 = new TCanvas("cComb", "Combined truth spectrum", 800, 600);
  c1->SetMargin(0.12, 0.05, 0.12, 0.08);
  c1->SetLogy();
  hComb->SetLineColor(kBlack);
  hComb->SetLineWidth(2);
  hComb->Draw("HIST");
  hComb->GetXaxis()->SetTitleSize(0.045);
  hComb->GetYaxis()->SetTitleSize(0.045);
  c1->SaveAs("LaBr_TruthSpectrum_Combined.pdf");

  // ---- Plot 2: photon gun only ---------------------------------------------
  TCanvas* c2 = new TCanvas("cGun", "Photon-gun truth spectrum", 800, 600);
  c2->SetMargin(0.12, 0.05, 0.12, 0.08);
  c2->SetLogy();
  hGun->SetLineColor(kBlue);
  hGun->SetLineWidth(2);
  hGun->Draw("HIST");
  hGun->GetXaxis()->SetTitleSize(0.045);
  hGun->GetYaxis()->SetTitleSize(0.045);
  {
    TLatex* l = new TLatex(); l->SetNDC(); l->SetTextSize(0.035);
    l->DrawLatex(0.50, 0.85, "Signal: 1809 keV #gamma");
    l->DrawLatex(0.50, 0.78, Form("Events: %ld", nGunEvts));
  }
  c2->SaveAs("LaBr_TruthSpectrum_PhotonGun.pdf");

  // ---- Plot 3: La-138 only --------------------------------------------------
  TCanvas* c3 = new TCanvas("cLa138", "La-138 truth spectrum", 800, 600);
  c3->SetMargin(0.12, 0.05, 0.12, 0.08);
  c3->SetLogy();
  hLa138->SetLineColor(kRed + 1);
  hLa138->SetLineWidth(2);
  hLa138->Draw("HIST");
  hLa138->GetXaxis()->SetTitleSize(0.045);
  hLa138->GetYaxis()->SetTitleSize(0.045);
  {
    TLatex* l = new TLatex(); l->SetNDC(); l->SetTextSize(0.035);
    l->DrawLatex(0.45, 0.85, "EC: 1436/1468 keV");
    l->DrawLatex(0.45, 0.78, "#beta^{-}: 789 keV + continuum");
    l->DrawLatex(0.45, 0.71, Form("Events: %ld", nLa138Evts));
  }
  c3->SaveAs("LaBr_TruthSpectrum_La138.pdf");

  std::cout << "\nDone. Output files:" << std::endl;
  std::cout << "  LaBr_TruthSpectrum_Combined.pdf" << std::endl;
  std::cout << "  LaBr_TruthSpectrum_PhotonGun.pdf" << std::endl;
  std::cout << "  LaBr_TruthSpectrum_La138.pdf" << std::endl;

  f->Close();
}
