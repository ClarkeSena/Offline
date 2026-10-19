// Dumps the STM truth pulses (STMMCPulse, from the LaBrTruthPulses ->
// TagPileupPulses chain) out of an art file into a flat CSV, and locates each
// pulse on the concatenated ADC .bin written by BinaryWaveformWriter.
//
// Reads the art file through ROOT's emulated streamers, so it needs no Offline
// dictionaries and runs under a bare ROOT anywhere, not just on the gpvm.
//
// Usage:
//   root -l -b -q 'dumpLaBrTruthPulses.C("sim.Poisson...art","truth.csv")'
//
// Sample-index conventions, validated against labr_sim_20kHz_ADC.bin (96.3% of
// isolated >100 keV pulses land exactly on peakSample, 100% within one sample):
//   - Event N's 628 samples occupy .bin indices [(N-1)*628, N*628); sample j of
//     event N is at absolute time (N-1)*1695 + j*2.7 ns. Note the truth timeline
//     advances 1695.0 ns per event while 628 samples span 1695.6 ns, so a global
//     time*fADC mapping drifts ~2 samples over a 10k-event run - hence per-event.
//   - STMMCPulse::time is the pulse LEADING EDGE, not its peak; the waveform
//     maximum arrives tPeak ~ 50 ns later (see below).
//   - ADC bins integrate charge over the bin, so a feature at time t lands in the
//     bin whose CENTRE is nearest t, i.e. floor(t/tADC - 0.5).
//   - Pulses are negative-going (PMT anode); peakADC is the |amplitude|.
//
// Original author: Sena Clarke

void dumpLaBrTruthPulses(const char* fin, const char* fout,
                         const char* processName = "PoissonPhotonGunToLaBrWaveform",
                         long long nEvents = 10000) {
  gErrorIgnoreLevel = kFatal;
  TFile* f = TFile::Open(fin);
  if (!f || f->IsZombie()) { printf("cannot open %s\n", fin); return; }
  TTree* t = (TTree*)f->Get("Events");
  TString br = TString::Format("mu2e::STMMCPulses_TagPileup__%s.obj", processName);

  auto grab = [&](const char* m) {
    Long64_t n = t->Draw(TString::Format("%s.%s", br.Data(), m), "", "goff");
    return std::vector<double>(t->GetV1(), t->GetV1() + n);
  };

  auto time      = grab("_time");
  auto energy    = grab("_energy");
  auto primaryId = grab("_primaryId");
  auto source    = grab("_source");
  auto eventId   = grab("_eventId");
  auto isPileup  = grab("_isPileup");
  auto nOverlap  = grab("_nOverlap");
  auto dtNearest = grab("_dtNearest");
  auto window    = grab("_window");

  const double micropulseTime = 1695.0;        // ns, LaBrPileup.micropulseTime
  const double tADC           = 2.7;           // ns, 1e3 / 370.370370 MSPS
  const int    nADCs          = 628;           // (int)(micropulseTime/tADC + 1)
  const double peakADCPerKeV  = 1.7543859649;  // LaBrPileup.peakADCPerKeV

  // Peak of the single-pe response (1-e^{-t/tau1})^p * e^{-t/tau2}: setting the
  // derivative to zero gives a closed form, tPeak = 49.95 ns = 18.5 samples.
  const double tau1 = 22.3543, tau2 = 28.0007, p = 6.6596;
  const double tPeak = -tau1 * log(tau1 / (p * tau2 + tau1));

  auto binOf = [&](double tLocal) { return (long long)floor(tLocal / tADC - 0.5); };

  const long long nSamplesTotal = nEvents * nADCs;
  long long nCarry = 0, nOutside = 0, nPileup = 0;

  FILE* o = fopen(fout, "w");
  fprintf(o, "pulseIndex,eventId,primaryId,source,sourceLabel,"
             "edgeTime_ns,localTime_ns,edgeSample,peakTime_ns,peakSample,"
             "inBin,carryOver,energy_keV,peakADC,"
             "isPileup,nOverlap,dtNearest_ns,window_ns\n");
  for (size_t i = 0; i < time.size(); ++i) {
    long long ev   = (long long)llround(eventId[i]);
    double local   = time[i] - (ev - 1) * micropulseTime;
    long long ejs  = binOf(local);
    long long pjs  = binOf(local + tPeak);
    long long eidx = (ev - 1) * nADCs + ejs;   // absorbs carry-over automatically
    long long pidx = (ev - 1) * nADCs + pjs;
    int carry      = (pjs >= nADCs) ? 1 : 0;   // charge spills into the next micropulse
    int inBin      = (pidx >= 0 && pidx < nSamplesTotal) ? 1 : 0;
    nCarry += carry;
    nOutside += (1 - inBin);
    nPileup += (int)isPileup[i];
    fprintf(o, "%zu,%lld,%lld,%d,%s,%.4f,%.4f,%lld,%.4f,%lld,%d,%d,%.4f,%.3f,%d,%d,%.4f,%.4f\n",
            i, ev, (long long)llround(primaryId[i]), (int)source[i],
            ((int)source[i] == 0 ? "external" : ((int)source[i] == 1 ? "La138" : ((int)source[i] == 2 ? "Ac227" : "unknown"))),
            time[i], local, eidx, time[i] + tPeak, pidx, inBin, carry,
            energy[i], energy[i] * peakADCPerKeV,
            (int)isPileup[i], (int)llround(nOverlap[i]), dtNearest[i], window[i]);
  }
  fclose(o);
  printf("\ndumpLaBrTruthPulses: %zu pulses -> %s (pileup %lld, carryOver %lld, outsideBin %lld)\n",
         time.size(), fout, nPileup, nCarry, nOutside);
}
