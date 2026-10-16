// Writes a LaBr pipeline job fcl whose photon spectrum is the line mix in a histogram with one
// bin per line, labelled by energy in keV (e.g. stmdet/hLinesCut: "66","347","844","1809").
// Bin contents become the relative lineIntensities (only the ratios matter).
//   rateHz = 0 : one photon per event at t = 0 (clean, evenly spaced pulses)
//   rateHz > 0 : Poisson pulse train at that photon rate, uniform times in each micropulse
// The job writes the ADC waveform both to DigiLaBr/ttree (nts file) and, via BinaryWriter, to a
// headerless little-endian int16 stream labr_SSCLines_<tag>_ADC.bin (np.fromfile(f,"<i2")).
//
// Usage: root -l -b -q 'Offline/STMMC/src/makeLinesFcl.C("gun_check_FOILOFF_w.root","stmdet/hLinesCut","FOILOFF")'
//        root -l -b -q 'Offline/STMMC/src/makeLinesFcl.C("gun_check_FOILOFF_w.root","stmdet/hLinesCut","FOILOFF_20kHz",20000)'
//   -> LaBrSSCLines_<tag>.fcl in the current directory; run it with
//      mu2e -c $PWD/LaBrSSCLines_<tag>.fcl
// Original author: Sena Clarke

#include "TFile.h"
#include "TH1.h"
#include "TString.h"

#include <cstdio>
#include <cstdlib>
#include <ctime>

void makeLinesFcl(const char* inFile, const char* histName = "stmdet/hLinesCut", const char* tag = "FOILOFF", double rateHz = 0.0) {
  TFile* f = TFile::Open(inFile);
  if (!f || f->IsZombie()) { printf("cannot open %s\n", inFile); return; }
  TH1* h = dynamic_cast<TH1*>(f->Get(histName));
  if (!h) { printf("no TH1 %s in %s\n", histName, inFile); return; }

  TString E, I;
  double total = 0;
  for (int i = 1; i <= h->GetNbinsX(); ++i) total += h->GetBinContent(i);
  printf("%s %s:\n", inFile, histName);
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    const char* label = h->GetXaxis()->GetBinLabel(i);
    const double keV = atof(label);
    if (keV <= 0) { printf("bin %d label '%s' is not an energy in keV, stopping\n", i, label); return; }
    E += Form("%s%.6g", i > 1 ? ", " : "", keV / 1000.0);
    I += Form("%s%.6g", i > 1 ? ", " : "", h->GetBinContent(i));
    printf("  %7.1f keV  weight %10.4g  (%.1f%%)\n", keV, h->GetBinContent(i), 100 * h->GetBinContent(i) / total);
  }

  const TString emission = rateHz > 0
    ? TString(Form("LaBrRun.primary.emission : @local::LaBrPipeline.Emission.poisson\n"
                   "LaBrRun.primary.emission.rateHz : %.6g\n", rateHz))
    : TString("LaBrRun.primary.emission : @local::LaBrPipeline.Emission.single\n");
  printf("emission: %s\n", rateHz > 0 ? Form("Poisson, %g Hz", rateHz) : "one photon per event at t = 0");

  // art forbids '_' in process_name; keep the tag as-is for file names
  TString procTag(tag);
  procTag.ReplaceAll("_", "");

  char date[16];
  time_t now = time(nullptr);
  strftime(date, sizeof(date), "%Y-%m-%d", localtime(&now));
  TString out = Form("LaBrSSCLines_%s.fcl", tag);
  FILE* o = fopen(out.Data(), "w");
  fprintf(o,
    "# LaBr pipeline with the line mix from %s %s\n"
    "# (relative weights, HPGe+LaBr combined). Generated %s by makeLinesFcl.C.\n"
    "#include \"Offline/STMMC/fcl/LaBrPipelineProlog.fcl\"\n"
    "BEGIN_PROLOG\n"
    "LaBrRun.processName : \"LaBrSSCLines%s\"\n"
    "LaBrRun.nEvents : 100000\n"
    "LaBrRun.primary : @local::LaBrPipeline.Primary.spectrum\n"
    "LaBrRun.primary.spectrum.lineEnergies    : [ %s ]\n"
    "LaBrRun.primary.spectrum.lineIntensities : [ %s ]\n"
    "LaBrRun.primary.beam.radius : 3.5\n"
    "%s"
    "LaBrRun.sources : [ primary ]\n"
    "LaBrRun.endModules : [ TruthLaBr, BinaryWriter, fullOutput ]\n"
    "LaBrRun.artFile    : \"sim.LaBrSSCLines%s.art\"\n"
    "LaBrRun.histFile   : \"nts.LaBrSSCLines%s.root\"\n"
    "LaBrRun.binaryFile : \"labr_SSCLines_%s_ADC.bin\"\n"
    "END_PROLOG\n"
    "#include \"Offline/STMMC/fcl/LaBrPipelineJob.fcl\"\n",
    inFile, histName, date, procTag.Data(), E.Data(), I.Data(), emission.Data(), tag, tag, tag);
  fclose(o);
  printf("wrote %s\n", out.Data());
}
