#!/usr/bin/env python3
"""Cross-check of the Geant4 27Mg excess against other yields from the SAME capture-dump events.

If the 3.4x 27Mg excess were a wrong count of stops or captures, every per-stop / per-capture yield would be off
by the same factor. This script prints, for Al stopping-target stops, Geant4/measured for several yields at once.

Input: the creationLine capture-dump CSVs (see analyseCaptureDump.py). Rows grouped by (file, x, y, z) into stops;
DIO = any row with model "DIO"; capture = a nu_mu (pdg 14) from model "NuclearCapture"; residual = heaviest ion
among the capture secondaries.

usage: captureCrossCheck.py file1.csv [file2.csv ...]      (writes <file1>.crosscheck.json)

Measured values:
  captures/stop            0.609 +- 0.001  Suzuki, Measday & Roalsvig, PRC 35 (1987) 2212: tau(Al) = 864.0 ns vs
                                           tau_free = 2197 ns gives 1 - 864/2197; +-0.001 covers the tau error
  neutrons/capture (mean)  1.262 +- 0.059  MacDonald et al., Phys. Rev. 139 (1965) B1253
  Al 2p-1s X-ray/stop      0.798 +- 0.008  Measday et al., PRC 76 (2007) 035504 (Geant4 line 334.673 keV,
                                           model EMCascade)
  1808.7 keV gamma/capture 0.51  +- 0.05   Measday et al., PRC 76 (2007) 035504 (26Mg 2+ -> 0+)
  26Mg/capture (observed)  0.53  +- 0.05   Measday et al., PRC 76 (2007) 035504
  27Mg/capture             0.0990+- 0.0033 Mizuno et al., arXiv:2507.19753 (activation)
  844 keV/capture          0.0711+- 0.0024 = 0.0990 x 0.718 (ENSDF I_gamma of 27Mg beta decay); the Geant4
                                           value is 27Mg/capture x 0.718: the decay is not in this dump
Errors: binomial for fractions (captures/stop, 26Mg, 27Mg); std/sqrt(N) for mean multiplicities (neutrons, X-rays,
1809 gammas); ratio errors add both relative errors in quadrature.
"""
import csv, sys, json, collections, math

sys.dont_write_bytecode = True      # importing the sibling script must not leave __pycache__ in the repo
from analyseCaptureDump import MATERIAL, ion_name

E_XRAY = 334.673    # keV, Geant4 Al 2p-1s
E_1809 = 1808.73    # keV, 26Mg 2+ -> 0+
WIN_1809 = 1.0      # keV
# The 1809 keV gamma is emitted by a recoiling 26Mg and is Doppler broadened by up to ~10 keV, so the yield is
# also taken in a +-20 keV window minus the continuum in the two 20 keV sidebands either side (same total width).
WIN_WIDE = 20.0     # keV
IG_844 = 0.718      # 27Mg beta decay, I_gamma(843.76 keV)
MG26, MG27 = 1000120260, 1000120270


def frac(n, d):
    return (n / d, math.sqrt(n * max(d - n, 0) / d**3)) if d else (0., 0.)


def mean_err(xs):
    n = len(xs)
    if not n:
        return (0., 0.)
    m = sum(xs) / n
    return (m, math.sqrt(max(sum(x * x for x in xs) / n - m * m, 0.) / n))


def ratio(g, m):
    r = g[0] / m[0]
    return (r, r * math.hypot(g[1] / g[0] if g[0] else 0., m[1] / m[0]))


def main(files):
    stops = collections.defaultdict(list)
    for fn in files:
        with open(fn) as f:
            for r in csv.DictReader(f):
                if r["material"] != MATERIAL:
                    continue
                stops[(fn, r["x_mm"], r["y_mm"], r["z_mm"])].append(
                    (int(r["pdgId"]), r["model"], float(r["ekin_keV"])))

    nstop = len(stops)
    ndio = ncap = nneither = nnoion = nmg26 = nmg27 = 0
    xrays, neutrons, g1809, g1809w = [], [], [], []
    g1809_other = collections.Counter()     # 1809 gammas in capture stops from models other than NuclearCapture
    lines_mg26, lines_all = collections.Counter(), collections.Counter()
    mg26_with_1809 = 0
    for parts in stops.values():
        xrays.append(sum(1 for p, m, e in parts if p == 22 and m == "EMCascade" and abs(e - E_XRAY) < 0.2))
        if any(m == "DIO" for _, m, _ in parts):
            ndio += 1
            continue
        cap = [(p, e) for p, m, e in parts if m == "NuclearCapture"]
        if not any(p == 14 for p, _ in cap):
            nneither += 1
            continue
        ncap += 1
        ions = [p for p, _ in cap if p > 1000020040]
        heavy = max(ions, key=lambda p: (p // 10) % 1000) if ions else 0
        nnoion += not heavy
        nmg26 += heavy == MG26
        nmg27 += heavy == MG27
        neutrons.append(sum(1 for p, _ in cap if p == 2112))
        gam = [e for p, e in cap if p == 22]
        g1809.append(sum(1 for e in gam if abs(e - E_1809) < WIN_1809))
        g1809w.append(sum(1 for e in gam if abs(e - E_1809) < WIN_WIDE)
                      - sum(1 for e in gam if WIN_WIDE <= abs(e - E_1809) < 2 * WIN_WIDE))
        for p, m, e in parts:
            if p == 22 and m != "NuclearCapture" and abs(e - E_1809) < WIN_1809:
                g1809_other[m] += 1
        for e in gam:
            lines_all[round(e)] += 1
        if heavy == MG26:
            for e in gam:
                lines_mg26[round(e)] += 1
            mg26_with_1809 += any(abs(e - E_1809) < WIN_1809 for e in gam)

    mg27 = frac(nmg27, ncap)
    rows = [  # (key, label, geant4, measured)
        ("cap_per_stop", "captures / (captures + DIO)", frac(ncap, ncap + ndio), (0.609, 0.001)),
        ("n_per_cap", "neutrons / capture (mean)", mean_err(neutrons), (1.262, 0.059)),
        ("xray_per_stop", "Al 2p-1s X-rays (334.673) / stop", mean_err(xrays), (0.798, 0.008)),
        ("g1809_per_cap", "1808.7 keV gammas / capture (+-1 keV)", mean_err(g1809), (0.51, 0.05)),
        ("g1809w_per_cap", "1808.7 keV gammas / capture (+-20 keV - sidebands)", mean_err(g1809w), (0.51, 0.05)),
        ("mg26_per_cap", "26Mg / capture", frac(nmg26, ncap), (0.53, 0.05)),
        ("mg27_per_cap", "27Mg / capture", mg27, (0.0990, 0.0033)),
        ("g844_per_cap", "844 keV / capture (=27Mg x 0.718, decay not in dump)",
         (mg27[0] * IG_844, mg27[1] * IG_844), (0.0711, 0.0024)),
    ]
    table = {k: {"label": lab, "geant4": g, "measured": m, "ratio": ratio(g, m)} for k, lab, g, m in rows}

    print("== Al target (%s): stops %d, DIO %d, captures %d, neither %d, captures with no residual ion %d" % (
        MATERIAL, nstop, ndio, ncap, nneither, nnoion))
    print("%-54s %-18s %-18s %s" % ("quantity", "Geant4", "measured", "Geant4/measured"))
    for k, lab, g, m in rows:
        r = table[k]["ratio"]
        print("%-54s %.4f +- %.4f   %.4f +- %.4f   %.3f +- %.3f" % (lab, g[0], g[1], m[0], m[1], r[0], r[1]))
    print("1808.7 keV gammas in capture stops from other models:", dict(g1809_other) or "none")

    f26 = frac(mg26_with_1809, nmg26)
    print("\n26Mg captures (%d): fraction with a capture gamma within %.1f keV of %.2f: %.4f +- %.4f" % (
        nmg26, WIN_1809, E_1809, f26[0], f26[1]))
    top26 = [(k, v, v / nmg26) for k, v in lines_mg26.most_common(10)] if nmg26 else []
    topall = [(k, v, v / ncap) for k, v in lines_all.most_common(10)] if ncap else []
    print("  top capture-gamma energies in 26Mg captures (keV, n, per 26Mg capture):")
    for k, v, f in top26:
        print("    %6d  %6d  %.4f" % (k, v, f))
    print("  top capture-gamma energies, all captures (keV, n, per capture):")
    for k, v, f in topall:
        print("    %6d  %6d  %.4f" % (k, v, f))

    others = [table[k]["ratio"][0] for k, _, _, _ in rows if k not in ("mg27_per_cap", "g844_per_cap")]
    r27 = table["mg27_per_cap"]["ratio"]
    print("\nGeant4/measured: other yields %.2f - %.2f; 27Mg %.2f +- %.2f" % (min(others), max(others), r27[0], r27[1]))

    res = {"files": files, "material": MATERIAL, "stops": nstop, "DIO": ndio, "captures": ncap,
           "neither": nneither, "captures_no_residual_ion": nnoion, "table": table,
           "g1809_other_models": dict(g1809_other), "mg26_frac_with_1809": f26,
           "mg26_gamma_lines_top10": top26, "all_capture_gamma_lines_top10": topall,
           "ratio_range_excl_27Mg": (min(others), max(others))}
    out = files[0] + ".crosscheck.json"
    with open(out, "w") as f:
        json.dump(res, f, indent=1)
    print("wrote", out)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])
