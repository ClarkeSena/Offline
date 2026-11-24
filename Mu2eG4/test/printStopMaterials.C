// Print mu- stops per material from the Run1B target muon-gun test.
//   root -l -b -q 'Offline/Mu2eG4/test/printStopMaterials.C("nts.run1bTargetMuonGun.closed.root")'
#include "TFile.h"
#include "TH1.h"
#include <cstdio>

void printStopMaterials(const char* fname = "nts.run1bTargetMuonGun.closed.root") {
  TFile f(fname);
  TH1* h = nullptr;
  f.GetObject("AllMuonFinder/stopmat", h);
  if (!h) { printf("AllMuonFinder/stopmat not found in %s\n", fname); return; }
  double tot = h->GetEntries();
  printf("%s: %.0f mu- stops\n", fname, tot);
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    double n = h->GetBinContent(i);
    if (n > 0) printf("  %-28s %8.0f  %5.1f%%\n", h->GetXaxis()->GetBinLabel(i), n, 100. * n / tot);
  }
}
