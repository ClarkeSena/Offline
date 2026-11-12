#!/usr/bin/env python3
# Table figure of the first MDC2025ai VD resampler run (s1 30743117 + s2 30750292):
# particles generated at VD116, how many reached the spot-size collimator (VD101), and the
# STM crystal hits with the VD116 particle each one started from.
# Numbers are from VDResampler_RunSummary.md section 4.3 and the crystal dump breakdown.
# Usage: python3 plotFirstRunStats.py  ->  FirstRun_VDResampler_Stats.png / .pdf
# Original author: Sena Clarke

import matplotlib.pyplot as plt

INK = "#1f2328"
MUTED = "#57606a"
RULE = "#d0d7de"
HEAD_BG = "#f3f4f6"
ACCENT_BG = "#e8f0fb"   # highlights the three rows asked for
SURFACE = "#ffffff"

NGEN = 1_040_000

flow_cols = ["Stage", "Plane", "Events", "Hits", "Fraction of\ngenerated events"]
flow_rows = [
    ["Generated at VD116", "VD116  STM_UpStrLarge  (z = 37700 mm, r = 2 m)", "1,040,000", "-", "1"],
    ["Crossings of VD116 in Geant4", "VD116", "-", "1,091,495", "-"],
    ["Downstream of field-of-view collimator", "VD100", "-", "4,159", "-"],
    ["Reached the spot-size collimator (SSC)", "VD101  STM_SpotSizeCollUpStr", "2,547", "2,967",
     f"{2547 / NGEN:.2e}"],
    ["Downstream of SSC", "VD88", "-", "3", "-"],
    ["LaBr front face", "VD89", "3", "3", f"{3 / NGEN:.1e}"],
    ["HPGe front face", "VD90", "1", "1", f"{1 / NGEN:.1e}"],
    ["Energy deposited in an STM crystal", "STMDet (HPGe or LaBr)", "7", "75 steps",
     f"{7 / NGEN:.1e}"],
    ["1809 keV gamma (+-1 keV) at VD116", "VD116 (~132 continuum + 14+-12 line)", "146", "-", f"{146 / NGEN:.1e}"],
    ["1809 keV gamma (+-1 keV) at VD101", "VD101", "1", "1", f"{1 / NGEN:.1e}"],
    ["1809 keV gamma in the STM", "STMDet", "0", "0", "0"],
]
flow_highlight = {0, 3, 7, 8, 9, 10}

src_cols = ["Source model", "Share of VD116 hits", "Expected events of 1.04e6"]
src_rows = [
    ["EleBeam", "90.72 %", "~943,500"],
    ["Neutrals", "6.99 %", "~72,700"],
    ["MuBeam", "2.29 %", "~23,800"],
    ["TargetStops1809", "0.0045 %", "~47"],
]

hit_cols = ["#", "Particle at VD116", "Kinetic E", "Distance from\nSTM axis", "Route to the crystal"]
hit_rows = [
    ["1", "gamma", "46.5 MeV", "~35 mm", "VD101 -> VD88 -> VD90: straight down the HPGe line of sight"],
    ["2", "gamma", "5.3 MeV", "~40 mm", "VD89: LaBr front face"],
    ["3", "neutron", "62.2 MeV", "~380 mm", "VD89: LaBr front face; 79Br(n,2n)78Br in LaBr"],
    ["4", "gamma", "15.6 MeV", "~140 mm", "No front-face crossing: through shielding / side"],
    ["5", "gamma", "7.9 MeV", "~1.29 m", "No front-face crossing: through shielding / side"],
    ["6", "gamma", "5.6 MeV", "~0.95 m", "No front-face crossing: through shielding / side"],
    ["7", "gamma", "0.5 MeV", "~1.26 m", "No front-face crossing: scattered in via shielding / side"],
]


def draw_table(ax, cols, rows, widths, highlight=(), fontsize=9.5):
    ax.axis("off")
    t = ax.table(cellText=rows, colLabels=cols, colWidths=widths, loc="upper left",
                 cellLoc="left", colLoc="left")
    t.auto_set_font_size(False)
    t.set_fontsize(fontsize)
    t.scale(1, 1.55)
    for (r, c), cell in t.get_celld().items():
        cell.set_edgecolor(RULE)
        cell.set_linewidth(0.6)
        cell.PAD = 0.03
        if r == 0:
            cell.set_facecolor(HEAD_BG)
            cell.get_text().set_color(INK)
            cell.get_text().set_fontweight("bold")
            cell.set_height(cell.get_height() * 1.35)
        else:
            cell.set_facecolor(ACCENT_BG if (r - 1) in highlight else SURFACE)
            cell.get_text().set_color(INK)
            if (r - 1) in highlight and c in (0, 2):
                cell.get_text().set_fontweight("bold")
    return t


fig = plt.figure(figsize=(13, 10.4), facecolor=SURFACE)
gs = fig.add_gridspec(4, 1, height_ratios=[0.5, 3.4, 1.45, 2.6], hspace=0.18,
                      left=0.03, right=0.97, top=0.965, bottom=0.115)

ax0 = fig.add_subplot(gs[0]); ax0.axis("off")
ax0.text(0, 0.95, "First VD resampler run: VD116 -> SSC -> STM crystals", fontsize=16,
         fontweight="bold", color=INK, va="top")
ax0.text(0, 0.30,
         "MDC2025ai models.  s1 generation job 30743117 (20 x 50k + 2 x 20k test events), "
         "s2 Geant4 job 30750292.  1.04e6 events = 7.69e9 POT (7396 POT/event), "
         "about 0.3 ms of beam.",
         fontsize=10, color=MUTED, va="top")

ax1 = fig.add_subplot(gs[1])
ax1.set_title("Particle flow (measured)", loc="left", fontsize=12, fontweight="bold", color=INK, pad=4)
draw_table(ax1, flow_cols, flow_rows, [0.27, 0.33, 0.10, 0.11, 0.14], flow_highlight)

ax2 = fig.add_subplot(gs[2])
ax2.set_title("What was generated at VD116, by source (expected from the hit-summary weights, not counted)",
              loc="left", fontsize=12, fontweight="bold", color=INK, pad=4)
draw_table(ax2, src_cols, src_rows, [0.22, 0.18, 0.22])

ax3 = fig.add_subplot(gs[3])
ax3.set_title("The 7 STM crystal-hit events and the VD116 particle each started from",
              loc="left", fontsize=12, fontweight="bold", color=INK, pad=4)
draw_table(ax3, hit_cols, hit_rows, [0.04, 0.14, 0.10, 0.12, 0.55], highlight=())

fig.text(0.03, 0.085,
         "Origin at VD116: 6 gamma + 1 neutron, 0 electrons.   Particles depositing energy in the crystals "
         "(75 steps): 46 e-, 25 gamma, 3 n, 1 ion (78Br).",
         fontsize=10, color=INK)
fig.text(0.03, 0.068,
         "Rate: 7 / 7.69e9 POT = 9.1e-10 per POT (+-38 % stat).  Only 2 of 7 came down the collimator line of sight; "
         "per-species counts at VD116 and VD101 not yet measured.\n"
         "Open check: VD101 / VD88 each show two 46.5474 MeV gamma entries in different events (possible duplicate).\n"
         "1809 keV gamma in +-1/10/100 keV windows: VD116 146/1331/13074, VD101 1/5/44 (continuum dominated at both).",
         fontsize=9, color=MUTED, va="top")

for ext in ("png", "pdf"):
    fig.savefig(f"FirstRun_VDResampler_Stats.{ext}", dpi=170, facecolor=SURFACE)
print("wrote FirstRun_VDResampler_Stats.png / .pdf")
