#!/usr/bin/env python3
# Where does a Mu2eG4Cuts creationLine cut act? Reads the CSV written with its dumpFile option (one line per
# selected track: Mu2e position, time, energy, creator process/model, volume, copy number, material) and prints
# counts by volume, material and process. With matplotlib available it also saves z-x and x-y scatter plots.
#
# Usage: python3 plotCreationPositions.py alxray_positions.csv [mg27_positions.csv ...]
#        (one PNG per input: <input>_positions.png)

import csv
import sys
from collections import Counter


def table(title, counter, total, top=25):
    print(f"\n  {title}")
    for key, n in counter.most_common(top):
        print(f"    {n:8d}  {100.0 * n / total:6.2f}%  {key}")
    if len(counter) > top:
        rest = sum(n for _, n in counter.most_common()[top:])
        print(f"    {rest:8d}  {100.0 * rest / total:6.2f}%  ({len(counter) - top} more)")


def summarise(fname):
    with open(fname, newline="") as f:
        rows = list(csv.DictReader(f))
    total = len(rows)
    print(f"\n=== {fname}: {total} tracks")
    if total == 0:
        return rows
    table("by volume (copy number)", Counter(f"{r['volume']} ({r['copyNo']})" for r in rows), total)
    table("by volume name", Counter(r["volume"] for r in rows), total)
    table("by material", Counter(r["material"] for r in rows), total)
    table("by creator process / model", Counter(f"{r['process']} / {r['model']}" for r in rows), total)
    return rows


def plot(fname, rows):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  (matplotlib not available: tables only)")
        return
    x = [float(r["x_mm"]) for r in rows]
    y = [float(r["y_mm"]) for r in rows]
    z = [float(r["z_mm"]) for r in rows]
    materials = [r["material"] for r in rows]
    order = [m for m, _ in Counter(materials).most_common(8)]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    for m in order:
        idx = [i for i, mm in enumerate(materials) if mm == m]
        ax1.scatter([z[i] for i in idx], [x[i] for i in idx], s=3, label=f"{m} ({len(idx)})")
        ax2.scatter([x[i] for i in idx], [y[i] for i in idx], s=3, label=f"{m} ({len(idx)})")
    ax1.set_xlabel("z (mm, Mu2e)")
    ax1.set_ylabel("x (mm, Mu2e)")
    ax1.set_title("creation points: z vs x")
    ax2.set_xlabel("x (mm, Mu2e)")
    ax2.set_ylabel("y (mm, Mu2e)")
    ax2.set_title("creation points: x vs y")
    ax1.legend(markerscale=4, fontsize=8)
    ax2.legend(markerscale=4, fontsize=8)
    fig.suptitle(fname)
    fig.tight_layout()
    out = fname.rsplit(".", 1)[0] + "_positions.png"
    fig.savefig(out, dpi=120)
    print(f"  saved {out}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__ if __doc__ else "usage: plotCreationPositions.py file.csv [...]")
    for name in sys.argv[1:]:
        plot(name, summarise(name))
