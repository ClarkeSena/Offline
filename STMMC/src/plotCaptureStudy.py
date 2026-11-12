#!/usr/bin/env python3
"""Figure for the 27Mg / muon-capture model study: Geant4 11.3.p02 Bertini (ShieldingM default) and
G4MuMinusCapturePrecompound against measurements, for muon capture in the Al stopping target.

Simulation numbers are the printed output of analyseCaptureDump.py (2026-11-10, 2 x 20k MuBeamCat MDC2025ab events
per model: 8725 Bertini and 8842 Precompound captures). Data: Measday et al. PRC 76, 035504 (2007) Table XI;
MacDonald et al. Phys. Rev. 139, B1253 (1965) (UCRL-11243 Table IV, folded with 54.5% neutron efficiency);
AlCap, arXiv:2110.10228 Table 5.

usage: plotCaptureStudy.py [out.png]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BERT, PREC, INK, MUTED, SURF = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#fcfcfb"
NB, NP = 8725, 8842

plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.grid": True, "grid.color": "#e4e3df", "grid.linewidth": 0.6, "axes.axisbelow": True,
                     "figure.facecolor": SURF, "axes.facecolor": SURF, "legend.frameon": False})
fig, axs = plt.subplots(2, 2, figsize=(12.5, 9))


def binerr(p, n):
    return np.sqrt(np.array(p) * (1 - np.array(p)) / n)


def bars(ax, x, b, p):
    w = 0.3
    ax.bar(x - w / 2 - 0.01, b, w, color=BERT, label="Geant4 Bertini (default)", yerr=binerr(b, NB),
           error_kw={"ecolor": INK, "lw": 0.8})
    ax.bar(x + w / 2 + 0.01, p, w, color=PREC, label="Geant4 Precompound", yerr=binerr(p, NP),
           error_kw={"ecolor": INK, "lw": 0.8})


# (a) residual nucleus per capture
ax = axs[0, 0]
cats = ["$^{27}$Mg\n(0n)", "$^{26}$Mg\n(1n)", "$^{25}$Mg\n(2n)", "$^{24}$Mg\n(3n)", "Na\n(A=23-26)", "Ne\n(A=21-23)"]
b = [0.3332, 0.3109, 0.0612, 0.1080, 0.1108, 0.0421]
p = [0.0144, 0.5102, 0.1884, 0.0788, 0.0746, 0.0875]
obs, eobs = [0.10, 0.53, 0.07, 0.02, 0.02, 0.01], [0.01, 0.05, 0.01, 0, 0, 0]
tot = [0.13, 0.61, 0.12, 0.06, 0.05, 0.03]
x = np.arange(len(cats))
bars(ax, x, b, p)
ax.errorbar(x, obs, yerr=eobs, fmt="o", color=INK, mec=SURF, mew=1.2, ms=8, capsize=3, label="Measday 2007, observed", zorder=5)
ax.plot(x, tot, "D", mfc=SURF, mec=INK, ms=7, label="Measday 2007, total (estimate)", zorder=5)
ax.set_xticks(x, cats)
ax.set_ylabel("yield per muon capture")
ax.set_title("(a) Residual nucleus", loc="left", fontweight="bold")
ax.text(-0.45, 0.36, "Bertini 0.333\nmeasured 0.10", fontsize=9, color=INK, ha="left", va="bottom")
ax.set_ylim(0, 0.72)
ax.legend(loc="upper right")

# (b) neutron multiplicity folded with 54.5 % efficiency
ax = axs[0, 1]
x = np.arange(4)
b = [0.539, 0.313, 0.114, 0.034]
p = [0.378, 0.498, 0.109, 0.015]
d, ed = [0.449, 0.464, 0.052, 0.036], [0.027, 0.028, 0.013, 0.007]
bars(ax, x, b, p)
ax.errorbar(x, d, yerr=ed, fmt="o", color=INK, mec=SURF, mew=1.2, ms=8, capsize=3, label="MacDonald 1965", zorder=5)
ax.set_xticks(x, ["0", "1", "2", "3"])
ax.set_xlabel("detected neutrons (true multiplicity folded with 54.5% efficiency)")
ax.set_ylabel("fraction of captures")
ax.set_title("(b) Neutron multiplicity", loc="left", fontweight="bold")
ax.text(0.96, 0.62, "mean true multiplicity\nBertini 1.18 $\\pm$ 0.01\nPrecompound 1.39 $\\pm$ 0.01\n"
        "measured 1.26 $\\pm$ 0.06", transform=ax.transAxes, ha="right", va="top", fontsize=9, color=INK)
ax.legend(loc="upper right")

# (c) neutrino energy
ax = axs[1, 0]
hb = {0: 3, 2: 12, 4: 16, 6: 17, 8: 29, 10: 40, 12: 36, 14: 47, 16: 46, 18: 56, 20: 75, 22: 76, 24: 76, 26: 76,
      28: 82, 30: 71, 32: 94, 34: 94, 36: 93, 38: 97, 40: 114, 42: 90, 44: 95, 46: 92, 48: 82, 50: 112, 52: 105,
      54: 104, 56: 94, 58: 101, 60: 94, 62: 96, 64: 117, 66: 104, 68: 99, 70: 95, 72: 98, 74: 92, 76: 84, 78: 89,
      80: 65, 82: 85, 84: 98, 86: 129, 88: 121, 90: 386, 92: 853, 94: 1104, 96: 1261, 98: 1043, 100: 487}
hp = {32: 8, 34: 7, 36: 19, 38: 49, 40: 56, 42: 69, 44: 83, 46: 101, 48: 121, 50: 133, 52: 187, 54: 192, 56: 206,
      58: 258, 60: 256, 62: 256, 64: 307, 66: 331, 68: 339, 70: 376, 72: 398, 74: 416, 76: 398, 78: 436, 80: 446,
      82: 515, 84: 501, 86: 480, 88: 527, 90: 527, 92: 520, 94: 295, 96: 29}
edges = np.arange(0, 106, 2)
for h, n, c, lab in ((hb, NB, BERT, "Geant4 Bertini (default)"), (hp, NP, PREC, "Geant4 Precompound")):
    y = np.array([h.get(int(e), 0) / n for e in edges[:-1]])
    ax.stairs(y, edges, color=c, lw=2, label=lab)
ax.axvline(102.1 - 6.443, color=INK, lw=1, ls="--")
ax.annotate("right of the dashed line the nucleus\nis too cold to emit a neutron\n"
            "($E^* < S_n$ = 6.44 MeV): $^{27}$Mg", (102.1 - 6.443, 0.075), (48, 0.085), ha="center", va="center",
            fontsize=9, color=INK, arrowprops={"arrowstyle": "-", "color": MUTED, "lw": 0.8})
ax.set_xlabel("neutrino energy (MeV)   [nuclear excitation $E^* \\approx 102.1 - E_\\nu$]")
ax.set_ylabel("fraction of captures per 2 MeV")
ax.set_title("(c) Neutrino energy = how cold the nucleus is left", loc="left", fontweight="bold")
ax.legend(loc="upper left")

# (d) charged particles in the AlCap energy windows, simulation / data
ax = axs[1, 1]
labs = ["p\n3.5-20 MeV", "d\n4.5-17 MeV", "t\n5-17 MeV", "$\\alpha$\n15-20 MeV"]
data = np.array([0.02664, 0.00846, 0.00258, 0.00044])
edata = np.array([0.00082, 0.00026, 0.00010, 0.00009])
x = np.arange(4)
for vals, errs, c, lab, dx in (([0.0989, 0.0074, 0.0023, 0.0010], [0.0032, 0.0009, 0.0005, 0.0003], BERT,
                                "Geant4 Bertini (default)", -0.08),
                               ([0.0425, 0.0141, 0.0053, 0.0005], [0.0021, 0.0013, 0.0008, 0.0002], PREC,
                                "Geant4 Precompound", 0.08)):
    r = np.array(vals) / data
    er = r * np.sqrt((np.array(errs) / np.array(vals))**2 + (edata / data)**2)
    ax.errorbar(x + dx, r, yerr=er, fmt="o", color=c, ms=8, capsize=3, label=lab)
ax.axhline(1, color=INK, lw=1)
ax.text(-0.45, 0.93, "AlCap 2021", ha="left", va="top", fontsize=9, color=INK)
ax.set_xticks(x, labs)
ax.set_xlim(-0.5, 3.5)
ax.set_ylim(0, 4.5)
ax.set_ylabel("simulation / measurement")
ax.set_title("(d) Charged-particle emission per capture", loc="left", fontweight="bold")
ax.legend(loc="upper right")

fig.suptitle("Muon capture in the Al stopping target: Geant4 11.3 capture models against measurements",
             fontsize=13, fontweight="bold", x=0.01, ha="left")
fig.tight_layout(rect=(0, 0, 1, 0.97))
out = sys.argv[1] if len(sys.argv) > 1 else "CaptureStudy_Mg27.png"
fig.savefig(out, dpi=160)
print("wrote", out)
