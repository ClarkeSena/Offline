#!/usr/bin/env python3
# Species breakdown of the first MDC2025ai VD resampler run (s1 30743117, s2 30750292):
# what was generated on VD116, what crossed VD101 (spot-size collimator upstream face), and
# the 7 STM crystal events (their VD116 origin and the tracks that deposited energy).
# Numbers from species_breakdown.sh (vd116_species.txt, vd101_species.txt, stm_species.txt).
# Usage: python3 plotSpeciesBreakdown.py  ->  SpeciesBreakdown_FirstRun.png / .pdf
# Original author: Sena Clarke

import matplotlib.pyplot as plt
import numpy as np

INK, MUTED, GRID, SURFACE = "#1f2328", "#57606a", "#e6e8eb", "#ffffff"
# fixed categorical order, one colour per species on every panel
SPECIES = ["gamma", "e-", "n", "e+", "ion (78Br)"]
COLOR = {"gamma": "#2a78d6", "e-": "#eb6834", "n": "#1baf7a", "e+": "#eda100",
         "ion (78Br)": "#e87ba4"}

L1809 = "gamma at 1809 keV (+-1 keV)"
COLOR[L1809] = "#4a3aa7"
# All measured. VD116 from the full-precision GenParticleDump (vd116_dump.sh): +-1/10/100 keV =
# 146 / 1331 / 13074, i.e. ~132 continuum (sidebands) + 14 +- 12 line. VD101: 1 / 5 / 44.
C1809 = {"VD116\ngenerated": 146, "VD101 (SSC face)\ncrossings": 1}
EXPECTED_1809 = set()

STAGES = ["VD116\ngenerated", "VD101 (SSC face)\ncrossings", "STM events:\norigin at VD116",
          "STM: tracks\ndepositing energy"]
COUNTS = {
    "VD116\ngenerated": {"gamma": 770621, "e-": 222133, "n": 44909, "e+": 2337},
    "VD101 (SSC face)\ncrossings": {"gamma": 2536, "e-": 293, "n": 38, "e+": 10},
    "STM events:\norigin at VD116": {"gamma": 6, "n": 1},
    "STM: tracks\ndepositing energy": {"e-": 46, "gamma": 25, "n": 3, "ion (78Br)": 1},
}
TOTAL_LABEL = {
    "VD116\ngenerated": "1.04e6 particles",
    "VD101 (SSC face)\ncrossings": "2,877 in 2,547 events",
    "STM events:\norigin at VD116": "7 events",
    "STM: tracks\ndepositing energy": "75 tracks, 7 events",
}


def style(ax):
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=9.5, length=0)
    ax.set_axisbelow(True)


fig = plt.figure(figsize=(14, 9.2), facecolor=SURFACE)
gs = fig.add_gridspec(2, 2, width_ratios=[1.55, 1], height_ratios=[4.6, 1.5],
                      left=0.06, right=0.985, top=0.86, bottom=0.04, wspace=0.16, hspace=0.42)
fig.text(0.06, 0.965, "First VD resampler run: particle species from VD116 to the STM crystals",
         fontsize=16, fontweight="bold", color=INK, va="top")
fig.text(0.06, 0.925, "MDC2025ai models, 1.04e6 events = 7.69e9 POT.  s1 30743117 (generation), "
         "s2 30750292 (Geant4).  VD101 counts are first crossings per track, any direction.",
         fontsize=10, color=MUTED, va="top")

# ---- left: counts per species and stage (log scale, one axis)
ax = fig.add_subplot(gs[0, 0]); style(ax)
x = np.arange(len(STAGES)); w = 0.145
for j, sp in enumerate(SPECIES):
    for i, st in enumerate(STAGES):
        v = COUNTS[st].get(sp, 0)
        if v == 0:
            continue
        xi = x[i] + (j - 2.5) * w
        ax.bar(xi, v, width=w * 0.88, color=COLOR[sp], edgecolor=SURFACE, linewidth=1)
        ax.text(xi, v * 1.25, f"{v:,}", ha="center", va="bottom", fontsize=8, color=INK,
                rotation=90 if v >= 1000 else 0)
# 1809 keV photons (subset of gamma), last bar in each group; hatched = expected, not measured
for i, st in enumerate(STAGES):
    xi = x[i] + 2.5 * w
    v = C1809.get(st, 0)
    if v == 0:
        ax.text(xi, 0.9, "0", ha="center", va="bottom", fontsize=8, color=COLOR[L1809], fontweight="bold")
        continue
    exp = st in EXPECTED_1809
    ax.bar(xi, v, width=w * 0.88, color=SURFACE if exp else COLOR[L1809], edgecolor=COLOR[L1809],
           hatch="///" if exp else None, linewidth=1.2)
    ax.text(xi, v * 1.25, f"~{v} exp." if exp else f"{v}", ha="center", va="bottom", fontsize=8,
            color=INK, rotation=90 if exp else 0)
ax.set_yscale("log"); ax.set_ylim(0.8, 6e6)
ax.set_xticks(x); ax.set_xticklabels([f"{st}\n({TOTAL_LABEL[st].replace(chr(10), ' ')})" for st in STAGES],
                                     color=INK, fontsize=9.5)
ax.yaxis.grid(True, color=GRID, linewidth=0.8)
ax.set_ylabel("Count (log scale)", color=MUTED, fontsize=10)
ax.set_title("Counts by species", loc="left", fontsize=12, fontweight="bold", color=INK)
handles = [plt.Rectangle((0, 0), 1, 1, color=COLOR[s]) for s in SPECIES]
handles.append(plt.Rectangle((0, 0), 1, 1, facecolor=COLOR[L1809]))
ax.legend(handles, SPECIES + ["gamma 1809 +-1 keV"], ncol=3, frameon=False,
          loc="upper right", fontsize=9, labelcolor=INK, bbox_to_anchor=(1.0, 1.12))

# ---- right: composition (each stage sums to 100 %)
ax2 = fig.add_subplot(gs[0, 1]); style(ax2)
bottom = np.zeros(len(STAGES))
for sp in SPECIES:
    frac = np.array([COUNTS[st].get(sp, 0) / sum(COUNTS[st].values()) for st in STAGES]) * 100
    ax2.bar(x, frac, bottom=bottom, width=0.6, color=COLOR[sp], edgecolor=SURFACE, linewidth=2)
    for i, f in enumerate(frac):
        if f >= 4:
            ax2.text(x[i], bottom[i] + f / 2, f"{f:.0f}%" if f >= 10 else f"{f:.1f}%",
                     ha="center", va="center", fontsize=8.5, color="#ffffff", fontweight="bold")
    bottom += frac
ax2.set_ylim(0, 100); ax2.set_yticks([0, 25, 50, 75, 100])
ax2.set_yticklabels(["0", "25", "50", "75", "100%"])
ax2.set_xticks(x); ax2.set_xticklabels(STAGES, color=INK, fontsize=8.5)
ax2.yaxis.grid(True, color=GRID, linewidth=0.8)
ax2.set_title("Composition at each stage", loc="left", fontsize=12, fontweight="bold", color=INK)

# ---- bottom: transmission per generated species
ax3 = fig.add_subplot(gs[1, :]); ax3.axis("off")
gen = COUNTS["VD116\ngenerated"]; v101 = COUNTS["VD101 (SSC face)\ncrossings"]
orig = COUNTS["STM events:\norigin at VD116"]
rows = []
for sp in ["gamma", "e-", "n", "e+"]:
    o = orig.get(sp, 0)
    stm = f"{o / gen[sp]:.1e}" if o else f"0  (< {2.3 / gen[sp]:.0e}, 90% CL)"
    rows.append([sp, f"{gen[sp]:,}", f"{v101.get(sp, 0):,}", f"{v101.get(sp, 0) / gen[sp]:.1e}",
                 str(o), stm])
rows.append(["gamma 1809\n+-1 keV", "146 (14+-12 line\n+ ~132 continuum)", "1", "6.8e-3", "0", "0  (< 0.016)"])
t = ax3.table(cellText=rows, colWidths=[0.07, 0.12, 0.11, 0.13, 0.12, 0.17], loc="upper left",
              colLabels=["Species", "Generated at VD116", "Crossing VD101", "VD101 / generated*",
                         "STM events (origin)", "STM events / generated"], cellLoc="left")
t.auto_set_font_size(False); t.set_fontsize(9.5); t.scale(1, 1.35)
t[(5, 0)].set_height(t[(5, 0)].get_height() * 1.8)
for c in range(6):
    t[(5, c)].set_height(t[(1, 0)].get_height() * 1.8)
for (r, c), cell in t.get_celld().items():
    cell.set_edgecolor("#d0d7de"); cell.set_linewidth(0.6); cell.PAD = 0.04
    cell.set_facecolor("#f3f4f6" if r == 0 else SURFACE)
    cell.get_text().set_color(INK)
    if r == 0:
        cell.get_text().set_fontweight("bold")
ax3.set_title("Transmission by species", loc="left", fontsize=12, fontweight="bold", color=INK, pad=2)
ax3.text(0.75, 0.95,
         "* species AT VD101, which can be a secondary of a\ndifferent generated particle "
         "(e.g. e- from a converted\ngamma).  VD101 forward-going (pz>0): 2,464.\n"
         "1809 row, +-1/10/100 keV: VD116 146/1331/13074,\nVD101 1/5/44. Continuum under the line from the\n"
         "+-10 keV sidebands: ~132 at VD116, ~0.2 at VD101.",
         transform=ax3.transAxes, fontsize=8.5, color=MUTED, va="top", ha="left")

for ext in ("png", "pdf"):
    fig.savefig(f"SpeciesBreakdown_FirstRun.{ext}", dpi=170, facecolor=SURFACE)
print("wrote SpeciesBreakdown_FirstRun.png / .pdf")
