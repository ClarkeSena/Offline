// Absolute rate accounting for the calibration-source run.
//
// Reads the truth-pulse tree written by TagPileupPulses (one entry per primary that
// deposited energy in the LaBr crystal, tagged by source) and turns it into counts per
// second, so the simulated external-source rate can be compared with a bench measurement
// and the intrinsic La-138 rate can be validated rather than scaled by hand.
//
// Every rate printed is divided by activityScale, so the numbers are PHYSICAL rates
// regardless of how much the statistics were boosted in the fcl.
//
// Sources: 0 = external calibration source, 1 = La-138, 2 = Ac-227 chain (only present when
// LaBrIntrinsicSource is used). Rates are also given in fixed energy windows chosen to
// isolate each feature, so they can be compared window by window with a bench spectrum:
//   Ba X-ray (20-60 keV), beta+789 (700-1100), 1473 sum peak (1380-1560),
//   alpha band (1560-2600, Ac-227 only; La-138 contributes almost nothing here).
// sourceActivity_Bq should be the DECAY-CORRECTED activity the fcl printed at startup.
//
// Usage:
//   root -l -b -q 'analyseCalibrationRates.C("nts.CalibrationSource.root",100,37000)'
//     arg1 file, arg2 activityScale (must match the fcl), arg3 source activity [Bq],
//     arg4 threshold [keV], arg5 nEvents fallback if it cannot be derived from the file.
//
// Original author: Sena Clarke

#include <cstdio>

void analyseCalibrationRates(const char* fname,
                             double activityScale = 100.0,
                             double sourceActivity_Bq = 37000.0,
                             double threshold_keV = 20.0,
                             long long nEventsFallback = 100000) {
  const double micropulseTime = 1695.0;  // ns, LaBrPileup.micropulseTime
  const int    nADCs          = 628;     // samples per micropulse at 370.37 MSPS

  TFile* f = TFile::Open(fname);
  if (!f || f->IsZombie()) { printf("cannot open %s\n", fname); return; }

  TTree* t = (TTree*)f->Get("TagPileup/pulses");
  if (!t) { printf("no TagPileup/pulses tree in %s (was makeTTree set?)\n", fname); return; }

  // Livetime. Derive the event count from the digitizer tree when it is present: it fills
  // exactly nADCs rows per micropulse, so this counts events that produced no pulse at all,
  // which the pulse tree cannot do.
  long long nEvents = nEventsFallback;
  TTree* dig = (TTree*)f->Get("DigiLaBr/ttree");
  if (dig && dig->GetEntries() % nADCs == 0) {
    nEvents = dig->GetEntries() / nADCs;
  } else {
    printf("NOTE: could not derive the event count from DigiLaBr/ttree; using %lld.\n"
           "      If that is wrong every rate below is wrong by the same factor.\n", nEvents);
  }
  const double livetime_s = nEvents * micropulseTime * 1e-9;

  double energy = 0, time = 0, dtNearest = 0, window = 0;
  int    source = 0, isPileup = 0, nOverlap = 0;
  unsigned eventId = 0, primaryId = 0;
  t->SetBranchAddress("energy",    &energy);
  t->SetBranchAddress("time",      &time);
  t->SetBranchAddress("source",    &source);
  t->SetBranchAddress("eventId",   &eventId);
  t->SetBranchAddress("primaryId", &primaryId);
  t->SetBranchAddress("isPileup",  &isPileup);
  t->SetBranchAddress("nOverlap",  &nOverlap);
  t->SetBranchAddress("dtNearest", &dtNearest);
  t->SetBranchAddress("window",    &window);

  // 8 keV bins, matching plotLaBrWaveform.C, extended to 2600 keV so the Co-60 sum peak
  // at 2505 keV is on scale.
  const int nBins = 375; const double eMin = 0.0, eMax = 3000.0;
  TH1D* hSrc = new TH1D("hSrc", "LaBr truth-pulse energy;Energy [keV];counts / s / 8 keV", nBins, eMin, eMax);
  TH1D* hLa  = new TH1D("hLa",  "La-138 intrinsic;Energy [keV];counts / s / 8 keV",        nBins, eMin, eMax);
  TH1D* hAc  = new TH1D("hAc",  "Ac-227 chain;Energy [keV];counts / s / 8 keV",            nBins, eMin, eMax);
  hSrc->SetDirectory(0); hLa->SetDirectory(0); hAc->SetDirectory(0);

  // Feature windows [keV] for the bench comparison.
  const int nWin = 4;
  const char* winName[nWin] = {"Ba X-ray 20-60", "beta+789 700-1100", "sum pk 1380-1560", "alpha 1560-2600"};
  const double winLo[nWin] = {20, 700, 1380, 1560};
  const double winHi[nWin] = {60, 1100, 1560, 2600};
  long long nWinCount[3][nWin] = {{0}};

  long long nSrc = 0, nLa = 0, nAc = 0, nUnk = 0, nSrcPile = 0, nLaPile = 0, nAcPile = 0;
  const Long64_t N = t->GetEntries();
  for (Long64_t i = 0; i < N; ++i) {
    t->GetEntry(i);
    if (energy < threshold_keV) continue;
    if (source == 0)      { ++nSrc; nSrcPile += isPileup; hSrc->Fill(energy); }
    else if (source == 1) { ++nLa;  nLaPile  += isPileup; hLa->Fill(energy); }
    else if (source == 2) { ++nAc;  nAcPile  += isPileup; hAc->Fill(energy); }
    else                  { ++nUnk; continue; }
    for (int w = 0; w < nWin; ++w)
      if (energy >= winLo[w] && energy < winHi[w]) ++nWinCount[source][w];
  }

  // Convert counts to physical rates: undo the statistics boost and divide by livetime.
  const double norm = 1.0 / (livetime_s * activityScale);
  const double rateSrc = nSrc * norm;
  const double rateLa  = nLa  * norm;
  const double rateAc  = nAc  * norm;
  hSrc->Scale(norm); hLa->Scale(norm); hAc->Scale(norm);

  printf("\n================ LaBr absolute rate accounting ================\n");
  printf("  file                      : %s\n", fname);
  printf("  events (micropulses)      : %lld\n", nEvents);
  printf("  livetime                  : %.6g s   (%lld x %.1f ns)\n", livetime_s, nEvents, micropulseTime);
  printf("  activityScale (divided out): %g\n", activityScale);
  printf("  energy threshold          : %.1f keV\n", threshold_keV);
  printf("  truth pulses in tree      : %lld\n", N);
  if (nUnk) printf("  WARNING: %lld pulses had source == unknown\n", nUnk);
  printf("\n  --- external calibration source (source == 0) ---\n");
  printf("    detected pulses         : %lld  (%lld tagged pileup, %.2f%%)\n",
         nSrc, nSrcPile, nSrc ? 100.0 * nSrcPile / nSrc : 0.0);
  printf("    DETECTED RATE           : %.4g counts/s   (above %.1f keV)\n", rateSrc, threshold_keV);
  if (sourceActivity_Bq > 0) {
    printf("    stated source activity  : %.4g Bq  (use the decay-corrected value)\n", sourceActivity_Bq);
    printf("    detected / decay        : %.4g   (total efficiency: solid angle x intrinsic\n",
           rateSrc / sourceActivity_Bq);
    printf("                                     x gamma yield per decay, all folded together)\n");
  }
  printf("\n  --- intrinsic La-138 (source == 1) ---\n");
  printf("    detected pulses         : %lld  (%lld tagged pileup, %.2f%%)\n",
         nLa, nLaPile, nLa ? 100.0 * nLaPile / nLa : 0.0);
  printf("    DETECTED RATE           : %.4g counts/s   (above %.1f keV)\n", rateLa, threshold_keV);
  printf("    The generator prints the La-138 activity it used at startup (~505 Bq computed,\n");
  printf("    ~495 Bq at the 1.425 Bq/cm^3 literature value, for a 347 cm^3 crystal). The\n");
  printf("    source is INSIDE the crystal, so the detected rate should be a large fraction.\n");
  if (nAc) {
    printf("\n  --- intrinsic Ac-227 chain (source == 2) ---\n");
    printf("    detected pulses         : %lld  (%lld tagged pileup, %.2f%%)\n",
           nAc, nAcPile, 100.0 * nAcPile / nAc);
    printf("    DETECTED RATE           : %.4g counts/s   (above %.1f keV)\n", rateAc, threshold_keV);
    if (rateLa > 0) printf("    La-138 / Ac-227 rate    : %.4g\n", rateLa / rateAc);
  }
  printf("\n  --- rate per feature window [counts/s]: compare with the bench spectrum ---\n");
  printf("    %-20s %12s %12s %12s %12s\n", "window [keV]", "source", "La-138", "Ac-227", "total");
  for (int w = 0; w < nWin; ++w) {
    const double r0 = nWinCount[0][w] * norm, r1 = nWinCount[1][w] * norm, r2 = nWinCount[2][w] * norm;
    printf("    %-20s %12.4g %12.4g %12.4g %12.4g\n", winName[w], r0, r1, r2, r0 + r1 + r2);
  }
  printf("    (on the bench, subtract a no-source run to remove residual room background)\n");
  printf("\n  --- the transfer-standard number ---\n");
  if (rateLa > 0) {
    printf("    source / La-138 rate    : %.4g\n", rateSrc / rateLa);
    printf("    Measure this same ratio on the bench with the real source at the same\n");
    printf("    distance. Agreement validates the La-138 activity; disagreement localises\n");
    printf("    the error, since the source side is known independently.\n");
  } else {
    printf("    no La-138 pulses (activityScale 0?), ratio not computed\n");
  }
  printf("===============================================================\n\n");

  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas("c", "calibration rates", 1000, 700);
  c->SetLogy();
  hSrc->SetLineColor(kAzure + 2); hSrc->SetLineWidth(2);
  hLa->SetLineColor(kOrange + 7); hLa->SetLineWidth(2);
  hAc->SetLineColor(kGreen + 2); hAc->SetLineWidth(2);
  double top = TMath::Max(TMath::Max(hSrc->GetMaximum(), hLa->GetMaximum()), hAc->GetMaximum());
  hSrc->SetMaximum(top * 5);
  hSrc->SetMinimum(TMath::Max(1e-6, norm * 0.5));
  hSrc->Draw("hist");
  hLa->Draw("hist same");
  if (nAc) hAc->Draw("hist same");
  TLegend* leg = new TLegend(0.55, 0.75, 0.88, 0.88);
  leg->AddEntry(hSrc, Form("calibration source (%.3g /s)", rateSrc), "l");
  leg->AddEntry(hLa,  Form("La-138 intrinsic (%.3g /s)", rateLa), "l");
  if (nAc) leg->AddEntry(hAc, Form("Ac-227 chain (%.3g /s)", rateAc), "l");
  leg->SetBorderSize(0);
  leg->Draw();
  c->SaveAs("CalibrationRates.pdf");
  printf("wrote CalibrationRates.pdf  (y axis is counts/s per 8 keV bin, physical rate)\n");

  f->Close();
}
