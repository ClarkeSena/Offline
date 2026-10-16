// Overlay the digitized energy spectra of the HPGe and LaBr STM detectors from a single
// combined simulation (FullSimToHPGeAndLaBrWaveforms.fcl). Both digitizers write a
// per-ADC-bin TTree whose chargeCollected branch sums, over an event, to the total charge
// the detector collected; converting that per-event sum to keV gives the energy estimate
// (the same approach plotLaBrWaveform.C uses for LaBr alone).
//   HPGe : chargeCollected is in Ge electron-hole pairs.
//          E[keV] = |sum_pairs| * epsilonGe[eV] / 1000, epsilonGe = 2.96 eV/pair.
//   LaBr : chargeCollected is in photoelectrons.
//          E[keV] = |sum_pe| / (lightYield * QE) * 1000, lightYield*QE = 63000*0.35.
// This is the pre-ADC charge estimator, so it carries the detector smearing (HPGe Fano,
// LaBr photostatistics + PMT excess noise) but not ADC quantization or electronic noise.
// Output: HPGeLaBr_EnergySpectra.pdf (overlay) and HPGeLaBr_EnergySpectra_panels.pdf.
// Usage: root -l -b -q 'plotHPGeAndLaBrSpectra.C("nts.FullSimHPGeLaBr.root")'
// Original author: Sena Clarke

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"

#include <map>
#include <cmath>
#include <iostream>

namespace {
  // Build a per-event energy spectrum from a digitizer TTree by summing |chargeCollected|
  // over the event and converting charge -> keV with the supplied factor.
  TH1D* makeSpectrum(TFile* f, const char* treePath, const char* histName,
                     const char* title, double keVPerCharge, int color) {
    TTree* t = (TTree*)f->Get(treePath);
    if (!t) {
      std::cerr << "WARNING: cannot find " << treePath
                << " (was the corresponding digitizer run with makeTTree:true?)" << std::endl;
      return nullptr;
    }
    double chargeCollected = 0;
    UInt_t eventId = 0;
    t->SetBranchAddress("chargeCollected", &chargeCollected);
    t->SetBranchAddress("eventId", &eventId);

    std::map<UInt_t, double> chargePerEvent; // eventId -> total |charge| this event
    Long64_t nEntries = t->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
      t->GetEntry(i);
      chargePerEvent[eventId] += std::fabs(chargeCollected);
    }

    const int    nBins   = 250;
    const double histMin = 0.0;
    const double histMax = 2000.0; // keV - covers the 1809 keV signal photopeak
    TH1D* h = new TH1D(histName, title, nBins, histMin, histMax);
    for (const auto& kv : chargePerEvent) {
      double e = kv.second * keVPerCharge;
      if (e > 0) h->Fill(e);
    }
    h->SetLineColor(color);
    h->SetLineWidth(2);
    std::cout << "  " << treePath << ": " << chargePerEvent.size()
              << " events, " << (long)h->GetEntries() << " with energy > 0" << std::endl;
    return h;
  }
}

void plotHPGeAndLaBrSpectra(const char* filename = "nts.FullSimHPGeLaBr.root") {
  gStyle->SetOptStat(0);

  TFile* f = TFile::Open(filename);
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open " << filename << std::endl;
    return;
  }

  // Charge -> keV conversions (see header comment).
  const double epsilonGe        = 2.96;            // eV per e-h pair
  const double hpgeKeVPerCharge = epsilonGe / 1e3; // keV per e-h pair
  const double lightYieldQE     = 63000.0 * 0.35;  // pe per MeV
  const double labrKeVPerCharge = 1e3 / lightYieldQE; // keV per pe

  std::cout << "Building spectra from " << filename << ":" << std::endl;
  TH1D* hHPGe = makeSpectrum(f, "DigiHPGe/ttree", "hHPGe",
    "STM digitized energy spectra;Energy deposited [keV];Counts / bin",
    hpgeKeVPerCharge, kRed + 1);
  TH1D* hLaBr = makeSpectrum(f, "DigiLaBr/ttree", "hLaBr",
    "STM digitized energy spectra;Energy deposited [keV];Counts / bin",
    labrKeVPerCharge, kBlue + 1);

  if (!hHPGe && !hLaBr) {
    std::cerr << "No detector TTrees found - nothing to plot." << std::endl;
    f->Close();
    return;
  }

  // ---- Plot 1: overlay (both detectors on one canvas) -----------------------
  TCanvas* c1 = new TCanvas("cOverlay", "HPGe + LaBr energy spectra", 900, 650);
  c1->SetMargin(0.12, 0.05, 0.12, 0.08);
  c1->SetLogy();

  double ymax = 1.0;
  if (hHPGe) ymax = std::max(ymax, hHPGe->GetMaximum());
  if (hLaBr) ymax = std::max(ymax, hLaBr->GetMaximum());

  // Draw the available histogram(s); first one sets the frame.
  TH1D* first = hLaBr ? hLaBr : hHPGe;
  first->SetMaximum(ymax * 2.0);
  first->SetMinimum(0.5);
  first->Draw("HIST");
  first->GetXaxis()->SetTitleSize(0.045);
  first->GetYaxis()->SetTitleSize(0.045);
  if (hLaBr && hHPGe) hHPGe->Draw("HIST SAME");

  TLegend* leg = new TLegend(0.58, 0.74, 0.93, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  if (hLaBr) leg->AddEntry(hLaBr, "LaBr3:Ce (scintillator)", "l");
  if (hHPGe) leg->AddEntry(hHPGe, "HPGe (semiconductor)", "l");
  leg->Draw();

  TLatex* tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.032);
  tex->DrawLatex(0.15, 0.85, "1.809 MeV #gamma into each detector");
  tex->DrawLatex(0.15, 0.80, "+ intrinsic ^{138}La in LaBr");

  c1->SaveAs("HPGeLaBr_EnergySpectra.pdf");
  std::cout << "Saved HPGeLaBr_EnergySpectra.pdf" << std::endl;

  // ---- Plot 2: side-by-side panels ------------------------------------------
  TCanvas* c2 = new TCanvas("cPanels", "HPGe + LaBr energy spectra (panels)", 1300, 550);
  c2->Divide(2, 1, 0.01, 0.01);

  c2->cd(1);
  gPad->SetMargin(0.13, 0.04, 0.12, 0.09);
  gPad->SetLogy();
  if (hHPGe) {
    hHPGe->SetTitle("HPGe digitized energy spectrum;Energy deposited [keV];Counts / bin");
    hHPGe->Draw("HIST");
  }

  c2->cd(2);
  gPad->SetMargin(0.13, 0.04, 0.12, 0.09);
  gPad->SetLogy();
  if (hLaBr) {
    hLaBr->SetTitle("LaBr3:Ce digitized energy spectrum;Energy deposited [keV];Counts / bin");
    hLaBr->Draw("HIST");
  }

  c2->SaveAs("HPGeLaBr_EnergySpectra_panels.pdf");
  std::cout << "Saved HPGeLaBr_EnergySpectra_panels.pdf" << std::endl;

  std::cout << "\nDone. Output files:" << std::endl;
  std::cout << "  HPGeLaBr_EnergySpectra.pdf" << std::endl;
  std::cout << "  HPGeLaBr_EnergySpectra_panels.pdf" << std::endl;

  f->Close();
}
