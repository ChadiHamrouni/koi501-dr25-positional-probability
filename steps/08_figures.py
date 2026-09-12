"""Step 8 - the two figures that answer questions the text can only assert.

Left: the r - J distribution of every tested neighbour, showing that the cut
used to call a colour impossible sits in a tail rather than through the bulk.

Right: separation against catalogued light share for the flagged neighbours,
with the region in which a scene model can prefer one marked, and every object
in it labelled by disposition.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from koi501 import config
from koi501.report import read

INK = "#161922"
MUTED = "#697285"
ACCENT = "#33468C"
BAD = "#A8382B"
OK = "#26624B"


def main() -> bool:
    step2 = read("02_colour_pathology")
    step5 = read("05_configuration")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.3))

    hist = step2["colour_histogram"]
    edges, counts = hist["edges"], hist["counts"]
    centres = [(a + b) / 2 for a, b in zip(edges, edges[1:])]
    colours = [BAD if c < config.COLOUR_IMPOSSIBLE else MUTED for c in centres]
    ax1.bar(centres, counts, width=0.24, color=colours, linewidth=0)
    ax1.set_yscale("log")
    ax1.axvline(config.COLOUR_IMPOSSIBLE, color=BAD, lw=1.2, ls="--")
    ax1.annotate(f"$r-J < {config.COLOUR_IMPOSSIBLE}$\nno star has this",
                 xy=(config.COLOUR_IMPOSSIBLE, max(counts) * 0.3),
                 xytext=(-8, 0), textcoords="offset points",
                 ha="right", va="center", fontsize=9, color=BAD)
    ax1.set_xlabel("$r - J$ from the catalogue's own two columns")
    ax1.set_ylabel(f"neighbours ({step2['n_pairs_tested']:,} tested)")
    ax1.set_title("The cut sits in a tail, not through the bulk",
                  fontsize=10, color=INK)
    ax1.set_xlim(-6, 6)

    sep_max, share_min = config.CONFIG_SEP_MAX, config.CONFIG_SHARE_MIN * 100
    ax2.add_patch(plt.Rectangle((0, share_min), sep_max, 100 - share_min,
                                facecolor=ACCENT, alpha=0.07, zorder=0))
    ax2.axvline(sep_max, color=ACCENT, lw=0.9, ls=":", zorder=1)
    ax2.axhline(share_min, color=ACCENT, lw=0.9, ls=":", zorder=1)

    inside = {(r["koi"], r["kic"]) for r in step5["table"]}
    outside = [r for r in step5["all_pairs"]
               if (r["koi"], r["kic"]) not in inside]
    ax2.scatter([r["sep_as"] for r in outside],
                [r["share"] * 100 for r in outside],
                c=MUTED, s=13, alpha=0.5, linewidth=0, zorder=2,
                label=f"flagged neighbour ({len(outside)})")

    style = {"FALSE POSITIVE": (BAD, "o"), "CONFIRMED": (ACCENT, "s"),
             "CANDIDATE": (OK, "*")}
    seen = set()
    for row in step5["table"]:
        colour, marker = style.get(row["disposition"], (MUTED, "o"))
        target = row["koi"] == config.KOI
        label = None
        if row["disposition"] not in seen:
            seen.add(row["disposition"])
            label = f"in the box, {row['disposition'].lower()}"
        ax2.scatter(row["sep_as"], row["share"] * 100, c=colour, marker=marker,
                    s=300 if target else 60, zorder=4 if target else 3,
                    edgecolor=INK if target else "none",
                    linewidth=1.0 if target else 0, label=label)
    ax2.annotate(config.KOI, xy=(6.74, 45), xytext=(-14, 14),
                 textcoords="offset points", fontsize=9, color=OK,
                 ha="right", fontweight="bold")
    ax2.annotate("close enough, and bright enough,\n"
                 "for the positional probability to prefer it",
                 xy=(0.6, 95), fontsize=8, color=ACCENT, va="top")
    ax2.set_xlabel("separation from the host (arcsec)")
    ax2.set_ylabel("neighbour's share of the catalogued light (%)")
    ax2.set_title(f"Only {len(step5['table'])} of "
                  f"{len(step5['all_pairs'])} land where it matters",
                  fontsize=10, color=INK)
    ax2.set_xlim(0, config.BLEND_RADIUS)
    ax2.set_ylim(0, 100)
    ax2.legend(frameon=False, fontsize=8, loc="upper right")

    for ax in (ax1, ax2):
        ax.tick_params(labelsize=9, colors=MUTED)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(MUTED)

    fig.tight_layout()
    config.FIGURES.mkdir(parents=True, exist_ok=True)
    for suffix in ("pdf", "png"):
        fig.savefig(config.FIGURES / f"colour_pathology.{suffix}", dpi=300)
    plt.close(fig)

    print(f"  wrote {config.FIGURES / 'colour_pathology.pdf'} (and .png)")
    below = sum(c for e, c in zip(hist["edges"], hist["counts"])
                if e < config.COLOUR_IMPOSSIBLE)
    print(f"  {below} of {step2['n_pairs_tested']:,} tested neighbours fall "
          f"below the cut ({100 * below / step2['n_pairs_tested']:.2f}%)")
    return True


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
