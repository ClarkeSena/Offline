// Compare the simulated Ac-227 alpha band (Ac227AlphaBand.fcl) with the bench measurement.
//
// Input: TruthLaBr/decays (LaBrTruthSpectrum), one entry per decay with its light-equivalent
// deposit (alpha steps already scaled by alphaQuenchingFactor/Slope). source == 2 is Ac-227.
// Energies are smeared with the measured resolution (FWHM/E at 1 MeV, scaling as 1/sqrt(E)),
// then counted in the 2.4-3.8 MeV band with the same linear sideband subtraction as
// labrSelfActivity.py (2.0-2.4 and 3.8-4.2 MeV).
//
//   T     = nEvents * windowTime * activityScale        (equivalent real time [s])
//   R_sim = N_band / T                                   (alphas/s at the fcl activity)
//
// It also prints the four band peak positions, to compare with the bench values
// 2.76 (Ra-223 + Th-227), 3.11 (Bi-211), 3.22 (Rn-219), 3.56 (Po-215) MeV.
//
// Usage:
//   root -l -b -q 'ac227AlphaBand.C("nts.Ac227AlphaBand.root", 20000, 1e5)'
//   optional args: windowTime [ns], rMeas, dRMeas, fwhmAt1MeV
// Original author: Sena Clarke

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TSpectrum.h"
#include "TStyle.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

void ac227AlphaBand(const char* filename = "nts.Ac227AlphaBand.root",
                    double nEvents = 0, double activityScale = 0, double windowTime = 1695.0,
                    double rMeas = 17.18, double dRMeas = 0.35, double fwhmAt1MeV = 0.026) {
  gStyle->SetOptStat(0);
  if (nEvents <= 0 || activityScale <= 0) {
    std::cerr << "Pass nEvents and activityScale (as in the fcl)" << std::endl;
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

  TH1D* h = new TH1D("h", "Ac-227 chain, light-equivalent deposit per decay;E [keV];decays / 10 keV",
                     500, 0, 5000);
  TRandom3 rng(12345);
  for (Long64_t i = 0; i < t->GetEntries(); ++i) {
    t->GetEntry(i);
    if (source != 2 || energy <= 0) continue;
    const double sigma = fwhmAt1MeV / 2.355 * std::sqrt(1000.0 * energy);
    h->Fill(rng.Gaus(energy, sigma));
  }

  auto integral = [&](double lo, double hi) {
    return h->Integral(h->FindBin(lo + 1e-6), h->FindBin(hi - 1e-6));
  };
  const double lo = 2400, hi = 3800;
  const double nGross = integral(lo, hi);
  const double nBg = 0.5 * (integral(2000, 2400) + integral(3800, 4200)) * (hi - lo) / 400.0;
  const double nBand = nGross - nBg;
  const double T = nEvents * windowTime * 1e-9 * activityScale;
  const double rSim = nBand / T, dRSim = std::sqrt(nGross) / T;

  std::cout << "\nequivalent time T = " << T << " s\n"
            << "alpha band " << lo << "-" << hi << " keV: " << nGross << " gross, " << nBg
            << " sideband -> " << nBand << "\n"
            << "simulated rate : " << rSim << " +- " << dRSim << " alphas/s\n"
            << "measured rate  : " << rMeas << " +- " << dRMeas << " alphas/s"
            << "  -> meas/sim = " << rMeas / rSim << " +- "
            << rMeas / rSim * std::hypot(dRMeas / rMeas, dRSim / rSim) << "\n";

  TH1D* hBand = (TH1D*)h->Clone("hBand");
  hBand->GetXaxis()->SetRangeUser(lo, hi);
  TSpectrum s(8);
  const int n = s.Search(hBand, 2, "nodraw goff", 0.05);
  std::vector<double> peaks(s.GetPositionX(), s.GetPositionX() + n);
  std::sort(peaks.begin(), peaks.end());
  std::cout << "band peaks [keV]:";
  for (double p : peaks) std::cout << " " << p;
  std::cout << "\nbench peaks      : 2760 3110 3220 3560\n" << std::endl;

  TCanvas* c = new TCanvas("c", "c", 900, 600);
  c->SetLogy();
  h->SetLineColor(kAzure + 2);
  h->GetXaxis()->SetRangeUser(0, 4500);
  h->Draw("hist");
  for (double e : {2760.0, 3110.0, 3220.0, 3560.0}) {
    TLine* l = new TLine(e, 0.5, e, h->GetMaximum());
    l->SetLineStyle(2);
    l->SetLineColor(kOrange + 7);
    l->Draw();
  }
  c->SaveAs("ac227AlphaBand.pdf");
}
