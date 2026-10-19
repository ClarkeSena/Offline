// La-138 internal full-absorption efficiency from La138InternalEfficiency.fcl, and the
// La-138 activity it implies for the bench-measured EC sum peak rate.
//
// Input: TruthLaBr/decays (LaBrTruthSpectrum), one entry per decay with the summed Geant4
// deposit of all its products. source == 1 is La-138.
//
// The EC sum peak is counted the same way as in labrSelfActivity.py: the 1380-1560 keV
// window minus a linear background, here taken from 1300-1360 and 1580-1640 keV sidebands.
// The beta branch (789 keV gamma + beta, endpoint 1044 keV) cannot reach the window, so
// every count in it is an EC decay. Truth energies are smeared with the measured
// resolution (FWHM/E at 1 MeV, scaling as 1/sqrt(E)) before counting.
//
//   eps_int = N_peak / (nDecays * ecFraction)
//   R_sim   = N_peak / T,  T = nDecays / A_sim   (EC peak rate at the simulated activity)
//   A_meas  = R_meas / (ecFraction * eps_int)
//
// nDecays and A_sim come from the job log ("La-138 decays: N" at endJob, and
// "La-138: ... -> X Bq" at construction). nDecays is the true generated count, so decays
// that deposited nothing (absent from the tree) are correctly counted as misses.
//
// Usage:
//   root -l -b -q 'la138InternalEfficiency.C("nts.La138InternalEfficiency.root", 171000, 505.2)'
//   optional args: rMeas (counts/s), dRMeas, fwhmAt1MeV, ecFraction
// Original author: Sena Clarke

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLine.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TStyle.h"

#include <cmath>
#include <iostream>

void la138InternalEfficiency(const char* filename = "nts.La138InternalEfficiency.root",
                             double nDecays = 0, double aSim = 0,
                             double rMeas = 63.438, double dRMeas = 0.666,
                             double fwhmAt1MeV = 0.026, double ecFraction = 0.655) {
  gStyle->SetOptStat(0);
  if (nDecays <= 0 || aSim <= 0) {
    std::cerr << "Pass nDecays and A_sim [Bq] from the job log (grep \"La-138\")" << std::endl;
    return;
  }
  TFile* f = TFile::Open(filename);
  if (!f || f->IsZombie()) { std::cerr << "Cannot open " << filename << std::endl; return; }
  TTree* t = (TTree*)f->Get("TruthLaBr/decays");
  if (!t) { std::cerr << "No TruthLaBr/decays in " << filename << std::endl; return; }

  double energy = 0;
  int source = -1;
  t->SetBranchAddress("energy", &energy);
  t->SetBranchAddress("source", &source);

  TH1D* hTruth = new TH1D("hTruth", "La-138 per-decay deposit;E [keV];decays / 2 keV", 1000, 0, 2000);
  TH1D* hSmear = new TH1D("hSmear", "La-138 per-decay deposit;E [keV];decays / 2 keV", 1000, 0, 2000);
  TRandom3 rng(12345);
  long nLa = 0;
  for (Long64_t i = 0; i < t->GetEntries(); ++i) {
    t->GetEntry(i);
    if (source != 1 || energy <= 0) continue;
    ++nLa;
    hTruth->Fill(energy);
    const double sigma = fwhmAt1MeV / 2.355 * std::sqrt(1000.0 * energy);
    hSmear->Fill(rng.Gaus(energy, sigma));
  }

  auto integral = [&](TH1D* h, double lo, double hi) {
    return h->Integral(h->FindBin(lo + 1e-6), h->FindBin(hi - 1e-6));
  };
  const double winLo = 1380, winHi = 1560;
  const double nWin = integral(hSmear, winLo, winHi);
  const double side = (integral(hSmear, 1300, 1360) + integral(hSmear, 1580, 1640)) / 120.0; // per keV
  const double nBg = side * (winHi - winLo);
  const double nPeak = nWin - nBg;
  const double dPeak = std::sqrt(nWin + nBg * (winHi - winLo) / 120.0);

  const double eps = nPeak / (nDecays * ecFraction);
  const double dEps = dPeak / (nDecays * ecFraction);
  const double T = nDecays / aSim;                  // equivalent real time [s]
  const double rSim = nPeak / T, dRSim = dPeak / T;
  const double aMeas = rMeas / (ecFraction * eps);
  const double dAMeas = aMeas * std::hypot(dRMeas / rMeas, dEps / eps);
  const double vol = M_PI * 3.81 * 3.81 * 7.62;     // cm^3

  std::cout << "\nLa-138 decays generated     : " << nDecays
            << "  (in tree with deposit: " << nLa << ")\n"
            << "EC window " << winLo << "-" << winHi << " keV : " << nWin
            << " counts, sideband bg " << nBg << " -> peak " << nPeak << " +- " << dPeak << "\n"
            << "eps_int (full EC absorption): " << eps << " +- " << dEps << "\n"
            << "simulated R_EC at " << aSim << " Bq: " << rSim << " +- " << dRSim << " counts/s\n"
            << "measured  R_EC              : " << rMeas << " +- " << dRMeas << " counts/s"
            << "  -> meas/sim = " << rMeas / rSim << "\n"
            << "implied A(La-138)           : " << aMeas << " +- " << dAMeas << " Bq ("
            << aMeas / vol << " Bq/cm^3; simulation uses " << aSim / vol << ")\n" << std::endl;

  TCanvas* c = new TCanvas("c", "c", 900, 600);
  c->SetLogy();
  hSmear->SetLineColor(kAzure + 2);
  hTruth->SetLineColor(kGray + 1);
  hSmear->Draw("hist");
  hTruth->Draw("hist same");
  for (double e : {winLo, winHi}) {
    TLine* l = new TLine(e, 0.5, e, hSmear->GetMaximum());
    l->SetLineStyle(2);
    l->SetLineColor(kOrange + 7);
    l->Draw();
  }
  TLegend* leg = new TLegend(0.15, 0.75, 0.5, 0.88);
  leg->AddEntry(hTruth, "Geant4 truth", "l");
  leg->AddEntry(hSmear, Form("smeared, FWHM %.1f%% at 1 MeV", 100 * fwhmAt1MeV), "l");
  leg->SetBorderSize(0);
  leg->Draw();
  c->SaveAs("la138InternalEfficiency.pdf");
}
