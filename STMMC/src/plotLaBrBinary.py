#!/usr/bin/env python3
"""Plot / export the raw LaBr ADC stream written by BinaryWaveformWriter.

The .bin file is headerless little-endian int16, one sample per 2.7 ns (370.37 MHz), 628
samples per 1695 ns micropulse, all micropulses back to back. Pulses are negative-going;
with NoiseSD = 0 the baseline is exactly 0. Energy (peak calibration) = 0.57 keV x |ADC peak|.

Examples
  # whole stream: min-envelope overview, every pulse visible (--plain: trace only)
  python3 plotLaBrBinary.py labr_SSCLines_FOILOFF_20kHz_ADC.bin --plain
  # a 10 ms section
  python3 plotLaBrBinary.py FILE.bin --start-us 50000 --length-us 10000 --plain
  # zoom: 50 us starting at 1000 us
  python3 plotLaBrBinary.py FILE.bin --start-us 1000 --length-us 50
  # micropulses 20..22
  python3 plotLaBrBinary.py FILE.bin --event 20 --nevents 3
  # pulse list (threshold crossing) to CSV, and the zoom range to .npz / .csv
  python3 plotLaBrBinary.py FILE.bin --pulses pulses.csv
  python3 plotLaBrBinary.py FILE.bin --start-us 1000 --length-us 50 --export zoom.csv
  # per-pulse energy spectrum (peak calibration), optionally with the truth spectrum
  python3 plotLaBrBinary.py FILE.bin --spectrum spectrum.png --truth truth.csv
  # overlay MC truth (CSV from TagPileup/pulses, see --help-truth)
  python3 plotLaBrBinary.py FILE.bin --event 20 --nevents 3 --truth truth.csv

In your own code:
  from plotLaBrBinary import load
  adc = load("FILE.bin")          # numpy memmap, int16, no copy
  t_us = np.arange(len(adc)) * TADC_NS * 1e-3
"""
import argparse
import sys

import numpy as np

FADC_MHZ = 370.370370370
TADC_NS = 1e3 / FADC_MHZ          # 2.7 ns
SAMPLES_PER_EVENT = 628           # int(1695 / 2.7) + 1, as the digitizer writes
MICROPULSE_NS = 1695.0            # truth timeline advances by this per event
PEAK_KEV_PER_ADC = 0.57           # LaBrDigitization.peakKeVPerADC

TRUTH_HELP = """Truth CSV: export TagPileup/pulses from the nts file with ROOT:
  root -l -b -q -e 'TFile f("nts.LaBrSSCLinesFOILOFF_20kHz.root"); TTree* t=(TTree*)f.Get("TagPileup/pulses"); FILE* o=fopen("truth.csv","w"); fprintf(o,"eventId,time,energy,source,isPileup\\n"); UInt_t ev; double tm,e; int s, p; t->SetBranchAddress("eventId",&ev); t->SetBranchAddress("time",&tm); t->SetBranchAddress("energy",&e); t->SetBranchAddress("source",&s); t->SetBranchAddress("isPileup",&p); for(Long64_t i=0;i<t->GetEntries();i++){t->GetEntry(i); fprintf(o,"%u,%.3f,%.3f,%d,%d\\n",ev,tm,e,s,p);} fclose(o);'
Truth times are absolute (event origin = (eventId-1) x 1695 ns) and mark the pulse LEADING edge;
the waveform peak comes ~50 ns later. The mapping to a sample is done per event, because
628 x 2.7 ns = 1695.6 ns, not 1695 ns (a global time -> sample mapping drifts)."""


def load(path):
    """The whole stream as a read-only int16 memmap (no copy, fast for large files)."""
    adc = np.memmap(path, dtype="<i2", mode="r")
    if len(adc) % SAMPLES_PER_EVENT:
        print(f"warning: {len(adc)} samples is not a multiple of {SAMPLES_PER_EVENT}", file=sys.stderr)
    return adc


def find_pulses(adc, threshold=-5):
    """Threshold crossings: one row per pulse (start, peak sample, peak ADC, energy keV)."""
    below = adc < threshold
    edges = np.flatnonzero(np.diff(below.astype(np.int8)))
    starts = edges[below[edges + 1]] + 1 if len(edges) else np.array([], dtype=int)
    if below[0]:
        starts = np.r_[0, starts]
    ends = edges[~below[edges + 1]] + 1 if len(edges) else np.array([], dtype=int)
    if below[-1]:
        ends = np.r_[ends, len(adc)]
    rows = []
    for s, e in zip(starts, ends):
        seg = np.asarray(adc[s:e])
        k = s + int(np.argmin(seg))
        rows.append((s, k, int(adc[k]), PEAK_KEV_PER_ADC * abs(int(adc[k]))))
    return np.array(rows, dtype=[("start", "i8"), ("peak", "i8"), ("adc", "i4"), ("energy_keV", "f8")])


def truth_samples(truth):
    """Map truth leading-edge times to stream sample indices, per event."""
    ev0 = truth["eventId"].astype(np.int64) - 1          # art events start at 1
    local = truth["time"] - ev0 * MICROPULSE_NS
    return ev0 * SAMPLES_PER_EVENT + np.floor(local / TADC_NS).astype(np.int64)   # the digitizer starts a pulse in bin floor(t / tADC)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("file")
    ap.add_argument("--start-us", type=float, help="zoom start [us of stream time]")
    ap.add_argument("--length-us", type=float, default=20.0, help="zoom length [us], default 20")
    ap.add_argument("--event", type=int, help="first micropulse to show (0-based)")
    ap.add_argument("--nevents", type=int, default=1, help="number of micropulses with --event")
    ap.add_argument("--threshold", type=int, default=-5, help="pulse threshold [ADC], default -5")
    ap.add_argument("--pulses", metavar="CSV", help="write the pulse list for the whole stream")
    ap.add_argument("--spectrum", metavar="PNG", help="plot the per-pulse energy spectrum and exit")
    ap.add_argument("--emax", type=float, default=0, help="spectrum upper edge [keV], default: highest pulse")
    ap.add_argument("--export", metavar="FILE", help="write the shown range (.csv or .npz)")
    ap.add_argument("--truth", metavar="CSV", help="overlay truth pulses (eventId,time,energy,source,isPileup)")
    ap.add_argument("--help-truth", action="store_true", help="how to make the truth CSV")
    ap.add_argument("--out", metavar="PNG", help="output image (default: derived from the file name)")
    ap.add_argument("--plain", action="store_true", help="ADC trace only: no truth markers, energy axis or micropulse lines")
    ap.add_argument("--show", action="store_true", help="open an interactive window as well")
    args = ap.parse_args()
    if args.help_truth:
        print(TRUTH_HELP)
        return

    adc = load(args.file)
    n = len(adc)
    print(f"{args.file}: {n} samples = {n / SAMPLES_PER_EVENT:.0f} micropulses = {n * TADC_NS * 1e-6:.3f} ms")

    if args.pulses:
        p = find_pulses(adc, args.threshold)
        with open(args.pulses, "w") as f:
            f.write("start_sample,peak_sample,peak_time_us,micropulse,peak_adc,energy_keV\n")
            for r in p:
                f.write(f"{r['start']},{r['peak']},{r['peak'] * TADC_NS * 1e-3:.4f},"
                        f"{r['peak'] // SAMPLES_PER_EVENT},{r['adc']},{r['energy_keV']:.1f}\n")
        print(f"{len(p)} pulses below {args.threshold} ADC -> {args.pulses}")

    if args.spectrum:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        p = find_pulses(adc, args.threshold)
        e = p["energy_keV"]
        emax = args.emax or 1.05 * max(e.max() if len(e) else 0, 100)
        bins = np.arange(0, emax + 8, 8.0)                  # 8 keV bins
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.hist(e, bins=bins, histtype="step", lw=1.5, label=f"digitized, per pulse ({len(e)})")
        if args.truth:
            tr = np.genfromtxt(args.truth, delimiter=",", names=True)
            ax.hist(tr["energy"], bins=bins, histtype="step", lw=1.2, ls="--",
                    label=f"truth, per pulse ({len(tr)})")
        over = (e > bins[-1]).sum()
        ax.set_yscale("log")
        ax.set_xlabel("Energy [keV]  (digitized: 0.57 keV x |peak ADC|)")
        ax.set_ylabel("Pulses / 8 keV")
        ax.set_title(f"{args.file}: pulse energy spectrum" + (f"  ({over} above {bins[-1]:.0f} keV)" if over else ""))
        ax.legend()
        fig.tight_layout()
        fig.savefig(args.spectrum, dpi=150)
        print(f"{len(e)} pulses, max {e.max():.0f} keV -> {args.spectrum}")
        return

    # Range to draw
    zoom = args.start_us is not None or args.event is not None
    if args.event is not None:
        lo = args.event * SAMPLES_PER_EVENT
        hi = lo + args.nevents * SAMPLES_PER_EVENT
    elif args.start_us is not None:
        lo = int(args.start_us * 1e3 / TADC_NS)
        hi = lo + int(args.length_us * 1e3 / TADC_NS)
    else:
        lo, hi = 0, n
    lo, hi = max(0, lo), min(n, hi)
    if hi <= lo:
        sys.exit("empty range")

    if args.export:
        seg = np.asarray(adc[lo:hi])
        idx = np.arange(lo, hi)
        if args.export.endswith(".npz"):
            np.savez_compressed(args.export, sample=idx, time_us=idx * TADC_NS * 1e-3, adc=seg)
        else:
            np.savetxt(args.export, np.c_[idx, idx * TADC_NS * 1e-3, seg], delimiter=",",
                       header="sample,time_us,adc", comments="", fmt=["%d", "%.4f", "%d"])
        print(f"samples {lo}-{hi} -> {args.export}")

    import matplotlib
    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(14, 5))
    # Long ranges are drawn as a min-envelope (most negative sample per block of samples), so
    # every pulse stays visible as a spike; short ranges are drawn sample by sample.
    nseg = hi - lo
    raw = nseg <= 20000
    long_axis = nseg * TADC_NS > 2e6                  # > 2 ms: time axis in ms
    tscale = 1e-6 if long_axis else 1e-3              # ns -> ms or us
    if raw:
        t = np.arange(lo, hi) * TADC_NS * tscale
        ax.plot(t, adc[lo:hi], lw=0.8, color="tab:blue")
        if not args.plain:
            for k in range(lo // SAMPLES_PER_EVENT, hi // SAMPLES_PER_EVENT + 1):   # micropulse edges
                ax.axvline(k * SAMPLES_PER_EVENT * TADC_NS * tscale, color="0.85", lw=0.6, zorder=0)
        title = f"samples {lo}-{hi}"
    else:
        width = 4000
        block = max(1, nseg // width)
        m = lo + (nseg // block) * block
        env = np.asarray(adc[lo:m]).reshape(-1, block).min(axis=1)
        t = (lo + np.arange(len(env)) * block + block / 2) * TADC_NS * tscale
        ax.fill_between(t, env, 0, step="mid", color="tab:blue", lw=0.4)
        what = "whole stream" if (lo == 0 and hi == n) else f"{nseg * TADC_NS * 1e-6:.2f} ms from {lo * TADC_NS * 1e-6:.2f} ms"
        title = f"{what}, min over {block} samples ({block * TADC_NS * 1e-3:.2f} us) per point"
    ax.set_xlabel("Stream time [ms]" if long_axis else "Stream time [us]")
    ax.set_xlim(t[0], t[-1])
    ax.set_ylabel("ADC")
    ax.ticklabel_format(axis="x", useOffset=False, style="plain")
    if not args.plain:
        sec = ax.secondary_yaxis("right", functions=(lambda a: -a * PEAK_KEV_PER_ADC, lambda e: -e / PEAK_KEV_PER_ADC))
        sec.set_ylabel("Peak energy [keV]")

    if args.truth and not args.plain:
        tr = np.genfromtxt(args.truth, delimiter=",", names=True)
        s = truth_samples(tr)
        sel = (s >= lo) & (s < hi)
        x = s[sel] * TADC_NS * tscale
        pile = tr["isPileup"][sel] > 0
        y = -tr["energy"][sel] / PEAK_KEV_PER_ADC
        ax.scatter(x[~pile], y[~pile], marker="v", color="tab:green", s=25, label="truth (isolated)", zorder=3)
        ax.scatter(x[pile], y[pile], marker="v", color="tab:red", s=25, label="truth (pileup)", zorder=3)
        ax.legend(loc="lower right")
        print(f"{sel.sum()} truth pulses in range ({pile.sum()} pileup)")

    ax.set_title(f"{args.file}: {title}")
    fig.tight_layout()
    out = args.out or args.file.rsplit(".", 1)[0] + ("_zoom.png" if zoom else "_overview.png")
    fig.savefig(out, dpi=150)
    print(f"plot -> {out}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
