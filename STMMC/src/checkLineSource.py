#!/usr/bin/env python3
"""Summarise a creationLine dump of the 347 or 844 keV line source (gen347check.fcl / gen844check.fcl).

usage: checkLineSource.py dump.csv nEvents
Prints photons per event for each line energy, where they were created, and their time distribution.
Expected (one target stop per event): 347 -> 0.798/event at 346.828 keV; 844 -> 0.0711 x 0.609 = 0.0433/event at
843.760 keV and 0.0279 x 0.609 = 0.0170/event at 1014.520 keV (ratio 0.392); all in TargetFoil_NN, StoppingTarget_Al.
"""
import csv, sys, collections, math
PERIOD, NB = 1695., 10  # micropulse (ns) = foldPeriod in BeamTo2VD844.fcl; flatness bins

if len(sys.argv) < 3:
    sys.exit(__doc__)
rows = list(csv.DictReader(open(sys.argv[1])))
nev = float(sys.argv[2])
print("%d photons in %s, %d events" % (len(rows), sys.argv[1], nev))
byE = collections.Counter(round(float(r["ekin_keV"]), 3) for r in rows)
for e, n in sorted(byE.items()):
    print("  %.3f keV: %d = %.4f +- %.4f per event" % (e, n, n / nev, math.sqrt(n) / nev))
es = sorted(byE)
if len(es) == 2:
    a, b = byE[es[0]], byE[es[1]]
    print("  ratio %.3f / %.3f = %.3f +- %.3f" % (es[1], es[0], b / a, b / a * math.sqrt(1. / a + 1. / b)))
print("creator process:", dict(collections.Counter(r["process"] for r in rows)))
print("material:", dict(collections.Counter(r["material"] for r in rows)))
vols = collections.Counter(r["volume"] for r in rows)
foils = [v for v in vols if v.startswith("TargetFoil")]
print("volumes: %d distinct, %d are TargetFoil_NN holding %.4f of the photons" % (
    len(vols), len(foils), sum(vols[v] for v in foils) / max(len(rows), 1)))
print("  photons per foil (foil number: n):", {int(v.split("_")[1]): vols[v] for v in sorted(foils)})
z = [float(r["z_mm"]) for r in rows]
print("z range %.1f - %.1f mm" % (min(z), max(z)))
for e in es:
    t = sorted(float(r["t_ns"]) for r in rows if round(float(r["ekin_keV"]), 3) == e)
    q = lambda f: t[min(int(f * len(t)), len(t) - 1)]
    print("time, %.3f keV: min %.1f, 10%% %.1f, median %.1f, 90%% %.1f, max %.1f ns" % (
        e, t[0], q(0.1), q(0.5), q(0.9), t[-1]))
    h = collections.Counter(int(x // 200) * 200 for x in t)
    print("  per 200 ns:", {k: h[k] for k in sorted(h)})
    # Folded lines (844, 1014) must lie in [0, PERIOD) and be flat; the 347 X-ray is not folded.
    nout = sum(1 for x in t if not 0. <= x < PERIOD)
    tin = [x for x in t if 0. <= x < PERIOD]
    hb = collections.Counter(int(x / (PERIOD / NB)) for x in tin)
    mu = len(tin) / float(NB)
    chi2 = sum((hb[i] - mu) ** 2 / mu for i in range(NB)) if mu > 0 else float("nan")
    print("  outside [0, %.0f) ns: %d of %d; inside, %d bins of %.1f ns: %s, chi2/ndf = %.1f/%d" % (
        PERIOD, nout, len(t), NB, PERIOD / NB, [hb[i] for i in range(NB)], chi2, NB - 1))
