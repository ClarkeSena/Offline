#!/usr/bin/env python3
"""Per-capture analysis of the Mu2eG4 creationLine capture dump (27Mg study).

Input: one or more CSVs written by the creationLine watch cut with process muMinusCaptureAtRest and all
secondary types (neutrinos, e-, gamma, nucleons, light ions, residual nuclei). All secondaries of one muon
stop share the same creation point, so rows are grouped by (x, y, z) into stops.

usage: analyseCaptureDump.py label out.json file1.csv [file2.csv ...]
"""
import csv, sys, json, collections, math

MATERIAL = "StoppingTarget_Al"
E0 = 102.1          # MeV, max neutrino energy for 27Al(mu-,nu)27Mg g.s. (m_mu - B_1s - nuclear mass difference)
SN_MG27 = 6.443     # MeV
LIGHT = {2112: "n", 2212: "p", 1000010020: "d", 1000010030: "t", 1000020030: "He3", 1000020040: "alpha"}
ELEM = {2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 10: "Ne", 11: "Na", 12: "Mg", 13: "Al"}


def ion_name(pdg):
    z = (pdg // 10000) % 1000
    a = (pdg // 10) % 1000
    return "%d%s" % (a, ELEM.get(z, "Z%d" % z))


def main():
    label, out = sys.argv[1], sys.argv[2]
    stops = collections.defaultdict(list)
    for fn in sys.argv[3:]:
        with open(fn) as f:
            for r in csv.DictReader(f):
                if r["material"] != MATERIAL:
                    continue
                stops[(fn, r["x_mm"], r["y_mm"], r["z_mm"])].append(
                    (int(r["pdgId"]), r["model"], float(r["ekin_keV"]) / 1000.0))

    nstop = len(stops)
    ndio = ncap = nneither = 0
    resid = collections.Counter()
    nmult = collections.Counter()
    charged = collections.Counter()        # any energy
    chargedWin = collections.Counter()     # AlCap energy windows
    win = {"p": (3.5, 20.), "d": (4.5, 17.), "t": (5., 17.), "alpha": (15., 20.)}
    enu, estar_mg27, nE = [], [], []
    enu_by_nmult = collections.defaultdict(list)
    gam_mg27 = collections.Counter()
    nxray = 0
    for parts in stops.values():
        models = set(m for _, m, _ in parts)
        nxray += sum(1 for p, m, e in parts if p == 22 and m == "EMCascade" and abs(e - 0.334673) < 2e-4)
        if "DIO" in models:
            ndio += 1
            continue
        cap = [(p, e) for p, m, e in parts if m == "NuclearCapture"]
        nus = [e for p, e in cap if p == 14]
        if not nus:
            nneither += 1
            continue
        ncap += 1
        ions = [p for p, e in cap if p > 1000020040]
        heavy = max(ions, key=lambda p: (p // 10) % 1000) if ions else 0
        name = ion_name(heavy) if heavy else "none"
        resid[name] += 1
        nn = sum(1 for p, e in cap if p == 2112)
        nmult[nn] += 1
        enu.append(nus[0])
        enu_by_nmult[min(nn, 4)].append(nus[0])
        nE += [e for p, e in cap if p == 2112]
        for p, e in cap:
            if p in LIGHT and p != 2112:
                k = LIGHT[p]
                charged[k] += 1
                if k in win and win[k][0] <= e <= win[k][1]:
                    chargedWin[k] += 1
        if name == "27Mg" and nn == 0:
            g = [e for p, e in cap if p == 22]
            estar_mg27.append(sum(g))
            for e in g:
                gam_mg27[round(e * 1000)] += 1

    def frac(n, d):
        return (n / d, math.sqrt(n * max(d - n, 0) / d**3)) if d else (0., 0.)

    res = {"label": label, "files": sys.argv[3:], "stops_Al_target": nstop, "DIO": ndio, "captures": ncap,
           "neither": nneither, "capture_fraction": frac(ncap, ncap + ndio),
           "xray_334p673_per_stop": frac(nxray, nstop)}
    res["residual_per_capture"] = {k: frac(v, ncap) + (v,) for k, v in resid.most_common()}
    res["neutron_multiplicity"] = {str(k): frac(v, ncap) for k, v in sorted(nmult.items())}
    mean_n = sum(k * v for k, v in nmult.items()) / ncap
    var_n = sum(k * k * v for k, v in nmult.items()) / ncap - mean_n**2
    res["mean_neutrons"] = (mean_n, math.sqrt(var_n / ncap))
    # MacDonald et al. Eq. 8: fold the true multiplicity with detection efficiency eps = 0.545
    eps = 0.545
    F = collections.defaultdict(float)
    for nu, v in nmult.items():
        for n in range(nu + 1):
            F[n] += v / ncap * math.comb(nu, n) * eps**n * (1 - eps)**(nu - n)
    res["folded_multiplicity_eps0.545"] = {str(k): F[k] for k in sorted(F)}
    res["charged_per_capture_any_energy"] = {k: frac(v, ncap) for k, v in charged.items()}
    res["charged_per_capture_AlCap_window"] = {k: frac(chargedWin[k], ncap) for k in win}
    estar = [E0 - e for e in enu]
    res["mean_Enu_MeV"] = sum(enu) / ncap
    res["mean_Estar_MeV(E0-Enu)"] = sum(estar) / ncap
    res["frac_Estar_below_Sn_6.443"] = frac(sum(1 for x in estar if x < SN_MG27), ncap)
    res["frac_Estar_below_1MeV"] = frac(sum(1 for x in estar if x < 1.0), ncap)
    res["frac_Enu_below_60MeV"] = frac(sum(1 for x in enu if x < 60.), ncap)
    res["mean_Enu_by_nmult"] = {str(k): (sum(v) / len(v), len(v)) for k, v in sorted(enu_by_nmult.items())}
    h = collections.Counter(int(e // 2) * 2 for e in enu)
    res["Enu_hist_2MeV"] = {str(k): h[k] for k in sorted(h)}
    hs = collections.Counter(min(int(x * 2) / 2, 30.0) for x in estar if x < SN_MG27 + 3)
    res["Estar_hist_0.5MeV_low"] = {str(k): hs[k] for k in sorted(hs)}
    hn = collections.Counter(min(int(e), 50) for e in nE)
    res["neutron_E_hist_1MeV"] = {str(k): hn[k] for k in sorted(hn)}
    res["neutron_frac_above_10MeV"] = frac(sum(1 for e in nE if e > 10.), len(nE))
    res["mean_neutron_E_MeV"] = sum(nE) / len(nE) if nE else 0.
    if estar_mg27:
        res["Mg27_sum_gamma_MeV_mean"] = sum(estar_mg27) / len(estar_mg27)
        res["Mg27_frac_no_gamma(gs)"] = frac(sum(1 for x in estar_mg27 if x == 0.), len(estar_mg27))
        hg = collections.Counter(int(x * 2) / 2 for x in estar_mg27)
        res["Mg27_sum_gamma_hist_0.5MeV"] = {str(k): hg[k] for k in sorted(hg)}
        res["Mg27_gamma_lines_keV_top"] = [(k, v, v / len(estar_mg27)) for k, v in gam_mg27.most_common(15)]
    json.dump(res, open(out, "w"), indent=1)

    print("== %s: Al target stops %d, DIO %d, captures %d, neither %d" % (label, nstop, ndio, ncap, nneither))
    print("capture fraction %.4f +- %.4f  (expected 0.609)" % res["capture_fraction"])
    print("334.673 keV X-rays per stop %.3f" % res["xray_334p673_per_stop"][0])
    print("residual nuclei per capture:")
    for k, v in list(res["residual_per_capture"].items())[:14]:
        print("   %-6s %.4f +- %.4f (%d)" % (k, v[0], v[1], v[2]))
    print("neutron multiplicity:", {k: "%.3f+-%.3f" % v for k, v in res["neutron_multiplicity"].items()})
    print("mean neutrons %.3f +- %.3f" % res["mean_neutrons"])
    print("folded (eps 0.545):", {k: round(v, 3) for k, v in res["folded_multiplicity_eps0.545"].items()})
    print("charged/capture any E:", {k: "%.4f" % v[0] for k, v in res["charged_per_capture_any_energy"].items()})
    print("charged/capture AlCap windows:", {k: "%.4f+-%.4f" % v for k, v in res["charged_per_capture_AlCap_window"].items()})
    print("mean Enu %.2f MeV, mean E* %.2f MeV" % (res["mean_Enu_MeV"], res["mean_Estar_MeV(E0-Enu)"]))
    print("frac E* < Sn(27Mg) %.4f +- %.4f; E* < 1 MeV %.4f; Enu < 60 MeV %.4f" % (
        res["frac_Estar_below_Sn_6.443"] + (res["frac_Estar_below_1MeV"][0], res["frac_Enu_below_60MeV"][0])))
    print("mean Enu by neutron multiplicity:", res["mean_Enu_by_nmult"])
    print("neutron mean E %.2f MeV, frac > 10 MeV %.3f" % (res["mean_neutron_E_MeV"], res["neutron_frac_above_10MeV"][0]))
    if estar_mg27:
        print("27Mg: mean sum-gamma %.2f MeV, no-gamma (g.s.) fraction %.3f" % (
            res["Mg27_sum_gamma_MeV_mean"], res["Mg27_frac_no_gamma(gs)"][0]))
        print("27Mg gamma lines (keV, n, per 27Mg):", [(k, v, round(f, 3)) for k, v, f in res["Mg27_gamma_lines_keV_top"]])
        print("27Mg sum-gamma hist:", res["Mg27_sum_gamma_hist_0.5MeV"])
    print("Enu hist (2 MeV bins):", res["Enu_hist_2MeV"])


if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    main()
