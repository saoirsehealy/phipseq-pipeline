#!/usr/bin/env python3
"""Genome / protein maps for a relabelled PhIP-Seq virus matrix.

One script in place of three. It replaces:

    dove/plot_genome_maps.py        % positive, whole + panel maps, per-ORF table
    plot_panel_genome_profiles.py   choice of y axis (pct / zscore / cpm / q75)
    dove/plot_genome_zscore.py      z on the y axis, with the mock-SD floor

and is virus-agnostic: nothing here is specific to SARS-CoV-2 or CCHFV.

WHAT GOES ON THE Y AXIS  (--measure, repeatable)

    pct      % of the group positive on that peptide. Breadth: how many people
             responded. Positive means z >= --z-threshold AND CPM >= --cpm-floor,
             which is the rule the ranking pipeline uses.
    zscore   z summarised across the group by --stat. Height: how strong the
             response was. A dashed line marks the threshold.
    cpm      raw CPM summarised by --stat, with a 25-75th percentile band.
    q75      75th percentile CPM minus the CPM floor, matching the pipeline's
             mock_q75 ranking mode.

THE MOCK SD FLOOR

    z is (CPM - mock mean) / mock SD. Many peptides have a mock SD of exactly
    zero -- 16% of the SARS-CoV-2 library -- which makes z undefined, and near-zero
    SDs send it to hundreds. The SD is therefore floored at the median of the
    non-zero mock SDs. The value used is printed and written into the figure.
    Pass --sd-floor 0 to disable, which reproduces a naive z exactly.

CHOOSING --stat

    mean     default; comparable between groups of different size
    median   conservative: most peptides sit near zero because only a minority
             of any group responds to any one peptide
    q90      what the strongest tenth reached; sensitive, but on a group of ~10
             it is effectively the maximum
    max      DO NOT use to compare groups. A maximum grows with the number of
             individuals and is set by single sticky peptides.

PEPTIDE NAMES

    <original>|<library tile>|<virus>|<ORF>|<genotype>|<start>-<end>

    The genotype may itself contain pipes (AmbiguousN:ref1|ref2|...), so the ORF
    is found by anchoring on the library tile and counting forward, never by
    counting back from the coordinates. Peptides carrying NOCOORD, and rows that
    were never relabelled, are reported and dropped rather than failing the run.

EXAMPLES

    # CCHFV M, one protein, amino-acid coordinates
    python3 plot_genome.py \
        --matrix cchfv_M_relabelled.csv \
        --panel  cchfv_m_relabelled_selected_panel.txt \
        --virus 'CCHFV M' --coordinate-level aa \
        --genome-length 1684 --genome-label NP_950235.1 \
        --mock-regex mock \
        --measure pct zscore \
        --group 'CCHFV=_CCHFV' --group 'Healthy=_neg'

    # SARS-CoV-2, ten ORFs, genome nucleotide coordinates
    python3 plot_genome.py \
        --matrix dove/prepared/cov2_relabelled.csv \
        --panel  dove/run_unsup_prepared/cov2_prep_selected_panel.txt \
        --virus SARS-CoV-2 --coordinate-level nt \
        --genome-length 29903 --genome-label NC_045512.2 \
        --mock-regex control_mockIP \
        --measure pct zscore \
        --group 'Vaccinated=^DOVE_' --group 'Infected=^ISARIC_' \
        --group 'Pre-pandemic=^Healthy_control_'
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
try:
    from plot_provenance import stamp
except Exception:                                    # provenance is optional
    stamp = None

INK, INK_2, GRID, SURFACE = "#0b0b0b", "#52514e", "#d8d7d2", "#fcfcfb"
ORF_COLOURS = ["#E8C33C", "#7B5EA7", "#7DC462", "#6E6E6E", "#E05C4B", "#1E7A46",
               "#5BC4C4", "#2B6CA3", "#E08A2E", "#B07AA1", "#9C6B4F", "#4A7DB5"]

TILE_RE = re.compile(r"^\d+-\d+(?:-\d+)?$")
COORD_RE = re.compile(r"^(\d+(?:\.\d+)?)-(\d+(?:\.\d+)?)$")


# ---------------------------------------------------------------- name parsing
def parse_names(names):
    """Return (orf, start, end) aligned to `names`; NaN/None where unusable."""
    orf, lo, hi = [], [], []
    for n in names.astype(str):
        parts = n.split("|")
        m = COORD_RE.match(parts[-1].strip()) if len(parts) >= 4 else None
        if not m:
            orf.append(None); lo.append(np.nan); hi.append(np.nan); continue

        tile = None
        for i in range(len(parts) - 2, -1, -1):
            if TILE_RE.match(parts[i].strip()):
                tile = i
                break
        if tile is None or tile + 2 > len(parts) - 2:
            orf.append(None); lo.append(np.nan); hi.append(np.nan); continue

        orf.append(parts[tile + 2].strip())
        lo.append(float(m.group(1))); hi.append(float(m.group(2)))
    idx = names.index
    return (pd.Series(orf, index=idx, dtype=object),
            pd.Series(lo, index=idx, dtype=float),
            pd.Series(hi, index=idx, dtype=float))


def group_spec(value):
    label, _, pattern = value.partition("=")
    label, pattern = label.strip(), pattern.strip()
    if not label or not pattern:
        raise argparse.ArgumentTypeError("--group must be 'Label=REGEX'")
    try:
        re.compile(pattern)
    except re.error as exc:
        raise argparse.ArgumentTypeError(f"bad regex {pattern!r}: {exc}")
    return label, pattern


def summarise(frame, stat):
    if stat == "mean":
        return frame.mean(axis=1)
    if stat == "median":
        return frame.median(axis=1)
    if stat == "q90":
        return frame.quantile(0.90, axis=1)
    return frame.max(axis=1)


# ---------------------------------------------------------------------- figure
def draw(measure, values, spread, orf, mid, groups, order, colour, args,
         out_dir, prefix, suffix, n_pep, sd_floor, n_floored, n_mock,
         panel_mid=None):
    """One figure: one stacked panel per group, bars coloured by ORF."""
    fig, axes = plt.subplots(len(groups), 1, sharex=True, sharey=True,
                             figsize=(13.2, 2.05 * len(groups) + 2.4),
                             squeeze=False)
    axes = axes[:, 0]
    fig.patch.set_facecolor(SURFACE)

    width = args.bar_width or (90.0 if args.coordinate_level == "nt" else 1.0)
    top = 0.0

    for ax, (label, cols) in zip(axes, groups):
        ax.set_facecolor(SURFACE)
        val = values[label]
        top = max(top, float(np.nanmax(val)) if len(val) else 0.0)
        for o in order:
            m = (orf == o).values
            if not m.any():
                continue
            ax.bar(mid[m], val[m], width=width, color=colour[o], linewidth=0,
                   zorder=3, label=o if ax is axes[0] else None)
            if spread is not None:
                lo_s, hi_s = spread[label]
                ax.vlines(mid[m], lo_s[m], hi_s[m], color=colour[o],
                          lw=0.8, alpha=0.55, zorder=2)
        # Panel peptides marked on the baseline, so the whole-library view shows
        # where the ranked panel sits relative to everything else.
        if panel_mid is not None and len(panel_mid):
            ax.plot(panel_mid, np.zeros(len(panel_mid)), marker="^",
                    linestyle="none", markersize=5, color=INK, zorder=6,
                    clip_on=False,
                    label="ranked panel peptide" if ax is axes[0] else None)

        if measure == "zscore":
            # z >= threshold is a PER-INDIVIDUAL rule. On a mean (or median) axis
            # it is not the level a real response has to reach: one positive donor
            # in a group of n contributes only 1/n of its z. Say so on the figure.
            ax.axhline(args.z_threshold, color=INK_2, lw=1.1, ls="--", zorder=4)
            if ax is axes[0]:
                note = (f"z = {args.z_threshold:g}"
                        if args.stat == "max"
                        else f"z = {args.z_threshold:g} per individual "
                             f"\u2014 not a threshold for the {args.stat}")
                ax.text(0.995, args.z_threshold, note, transform=ax.get_yaxis_transform(),
                        fontsize=8, color=INK_2, va="bottom", ha="right", zorder=5)
        ax.set_title(f"{label}   (n={len(cols)})", fontsize=11, color=INK, pad=6)
        ax.grid(axis="y", color=GRID, lw=0.6, zorder=0)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(GRID)
        ax.tick_params(colors=INK_2, labelsize=9)

    axes[0].set_ylim(0, 100 if measure == "pct" else max(top * 1.08, 1e-9))
    axes[0].legend(frameon=False, fontsize=8.5, ncol=min(len(order), 11),
                   loc="lower left", bbox_to_anchor=(0, 1.28),
                   columnspacing=1.0, handlelength=1.1)

    ylab = {
        "pct": f"Individuals positive (%)\n[z ≥ {args.z_threshold:g} "
               f"and CPM ≥ {args.cpm_floor:g}]",
        "zscore": f"{args.stat} z-score across the group\n"
                  f"[z = (CPM − mock mean) / mock SD]",
        "cpm": f"{args.stat} CPM across the group",
        "q75": f"75th-percentile CPM − {args.cpm_floor:g}",
    }[measure]
    # centred across all panels, so it stays on the figure whatever
    # the group count
    fig.supylabel(ylab, fontsize=10.5, color=INK_2, x=0.012)

    unit = "nt" if args.coordinate_level == "nt" else "aa"
    xlab = f"Genome location ({unit})"
    if args.genome_label:
        xlab += f"  ({args.genome_label})"
    axes[-1].set_xlabel(xlab, fontsize=10.5, color=INK_2)
    axes[-1].set_xlim(0, args.genome_length or float(np.nanmax(mid)) * 1.02)

    fig.suptitle(f"{args.virus} — {n_pep:,} peptides", fontsize=13.5,
                 color=INK, y=0.995)

    rect_bottom = 0.10
    if stamp is not None:
        fig.tight_layout(rect=[0, rect_bottom, 1, 0.90])
        note = {
            "pct": f"bar = % of that group with z ≥ {args.z_threshold:g} "
                   f"and CPM ≥ {args.cpm_floor:g}",
            "zscore": f"bar = {args.stat} z across the group; dashed line "
                      f"z = {args.z_threshold:g}",
            "cpm": f"bar = {args.stat} CPM; whisker = 25th–75th percentile",
            "q75": f"bar = 75th-percentile CPM minus {args.cpm_floor:g}",
        }[measure]
        if panel_mid is not None and len(panel_mid):
            note += f" | triangles = {len(panel_mid)} ranked panel peptides"
        stamp(fig, Path(args.matrix).name,
              f"z vs {n_mock} mock wells, SD floored at {sd_floor:.3g} "
              f"({n_floored:,} peptides at the floor)",
              extra=note)
        fig.tight_layout(rect=[0, rect_bottom + 0.035, 1, 0.90])
    else:
        fig.tight_layout(rect=[0, 0.02, 1, 0.90])

    png = out_dir / f"{prefix}_genome_{measure}{suffix}.png"
    fig.savefig(png, dpi=args.dpi, facecolor=SURFACE)
    plt.close(fig)
    print(f"Wrote {png}")
    return png


# ------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(
        description="Genome/protein maps for a relabelled PhIP-Seq virus matrix.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix", required=True,
                    help="relabelled CPM matrix from prepare_virus_matrix_v2.py")
    ap.add_argument("--panel", help="selected panel txt; adds a panel-only figure")
    ap.add_argument("--virus", default="virus", help="name used in titles")
    ap.add_argument("--out-dir", default=".", help="directory for figures and tables")
    ap.add_argument("--prefix", help="output filename prefix; defaults to --virus")
    ap.add_argument("--measure", nargs="+", default=["pct"],
                    choices=["pct", "zscore", "cpm", "q75"],
                    help="what goes on the y axis; repeatable")
    ap.add_argument("--stat", default="mean",
                    choices=["mean", "median", "q90", "max"],
                    help="how zscore and cpm are summarised across a group")
    ap.add_argument("--group", action="append", type=group_spec, required=True,
                    metavar="LABEL=REGEX", help="sample group; repeatable")
    ap.add_argument("--mock-regex", default="mock",
                    help="pattern identifying beads-only/mock columns")
    ap.add_argument("--coordinate-level", choices=["nt", "aa"], default="aa")
    ap.add_argument("--genome-length", type=float, help="x-axis maximum")
    ap.add_argument("--genome-label", help="reference accession for the x label")
    ap.add_argument("--z-threshold", type=float, default=3.5)
    ap.add_argument("--cpm-floor", type=float, default=10.0)
    ap.add_argument("--sd-floor", type=float, default=None,
                    help="floor for the mock SD; default is the median non-zero "
                         "mock SD. Pass 0 for a naive z.")
    ap.add_argument("--no-spread", action="store_true",
                    help="omit the 25-75th percentile whisker on cpm")
    ap.add_argument("--bar-width", type=float)
    ap.add_argument("--orf-order", help="comma-separated ORF order for the legend")
    ap.add_argument("--dpi", type=int, default=200)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = args.prefix or re.sub(r"[^\w.-]+", "_", args.virus).strip("_")

    # ---- load and parse -----------------------------------------------------
    d = pd.read_csv(args.matrix, index_col=0)
    names = pd.Series(d.index.astype(str), index=d.index)
    orf, lo, hi = parse_names(names)

    usable = orf.notna() & lo.notna() & hi.notna()
    if not usable.any():
        raise SystemExit(
            "ERROR: no peptide carries usable coordinates. Check the matrix came "
            "from prepare_virus_matrix_v2.py.")
    if (~usable).any():
        print(f"NOTE: dropping {int((~usable).sum()):,} of {len(d):,} peptides "
              "with no usable coordinates (NOCOORD or never relabelled).")
    d, names = d[usable.values], names[usable]
    orf, lo, hi = orf[usable], lo[usable], hi[usable]
    mid = ((lo + hi) / 2.0).values

    # ---- mocks and z --------------------------------------------------------
    mock_cols = [c for c in d.columns if re.search(args.mock_regex, str(c))]
    if not mock_cols:
        raise SystemExit(f"ERROR: no column matched --mock-regex {args.mock_regex!r}")
    mu = d[mock_cols].mean(axis=1)
    sd = d[mock_cols].std(axis=1)
    nonzero = sd[sd > 1e-6]
    sd_floor = (args.sd_floor if args.sd_floor is not None
                else (float(nonzero.median()) if len(nonzero) else 1.0))
    n_floored = int((sd <= sd_floor).sum())
    denom = np.maximum(sd, sd_floor) if sd_floor > 0 else sd.replace(0.0, np.nan)
    Z = d.sub(mu, axis=0).div(denom, axis=0).fillna(0.0)

    # ---- groups -------------------------------------------------------------
    groups = []
    for label, pattern in args.group:
        pat = re.compile(pattern)
        cols = [c for c in d.columns
                if pat.search(str(c)) and c not in mock_cols]
        if cols:
            groups.append((label, cols))
        else:
            print(f"NOTE: group {label!r} (/{pattern}/) matched no non-mock column.")
    if not groups:
        raise SystemExit("ERROR: no --group matched any non-mock column.")

    print(f"Matrix: {args.matrix}")
    print(f"  {len(d):,} peptides | {len(mock_cols)} mock wells | stat: {args.stat}")
    print(f"  mock SD floored at {sd_floor:.4g} "
          f"({n_floored:,} peptides at or below it, {100*n_floored/len(d):.1f}%)")
    print("  " + " | ".join(f"{lab} n={len(c)}" for lab, c in groups))

    order = list(pd.unique(pd.Series(orf).iloc[np.argsort(mid)]))
    if args.orf_order:
        wanted = [o.strip() for o in args.orf_order.split(",") if o.strip()]
        order = [o for o in wanted if o in order] + [o for o in order if o not in wanted]
    colour = {o: ORF_COLOURS[i % len(ORF_COLOURS)] for i, o in enumerate(order)}

    panel = set()
    if args.panel and Path(args.panel).exists():
        panel = {l.strip() for l in open(args.panel) if l.strip()}
        found = int(names.isin(panel).sum())
        print(f"  panel peptides found in matrix: {found:,} of {len(panel):,}")

    # ---- each measure, for the whole library and then the panel -------------
    views = [("", np.ones(len(d), bool))]
    if panel:
        keep = names.isin(panel).values
        if keep.any():
            views.append(("_panel", keep))
        else:
            print("  panel peptides not found in this matrix; skipping panel figure.")

    panel_positions = mid[names.isin(panel).values] if panel else np.array([])

    rows = []
    for suffix, mask in views:
        sub_orf, sub_mid = orf[mask], mid[mask]
        # only on the whole-library view; on the panel figure every bar is a
        # panel peptide, so the markers would say nothing
        marks = panel_positions if (suffix == "" and len(panel_positions)) else None
        for measure in args.measure:
            values, spread = {}, ({} if measure == "cpm" and not args.no_spread else None)
            for label, cols in groups:
                if measure == "pct":
                    pos = (Z[cols] >= args.z_threshold) & (d[cols] >= args.cpm_floor)
                    v = pos.mean(axis=1) * 100.0
                elif measure == "zscore":
                    v = summarise(Z[cols], args.stat)
                elif measure == "cpm":
                    v = summarise(d[cols], args.stat)
                    if spread is not None:
                        spread[label] = (d[cols].quantile(0.25, axis=1).values[mask],
                                         d[cols].quantile(0.75, axis=1).values[mask])
                else:
                    v = (d[cols].quantile(0.75, axis=1) - args.cpm_floor).clip(lower=0)
                values[label] = v.values[mask]

            draw(measure, values, spread, sub_orf, sub_mid, groups, order, colour,
                 args, out_dir, prefix, suffix, int(mask.sum()),
                 sd_floor, n_floored, len(mock_cols), panel_mid=marks)

            for label, _ in groups:
                v = values[label]
                for o in order:
                    m = (sub_orf == o).values
                    if not m.any():
                        continue
                    rows.append(dict(view="panel" if suffix else "all",
                                     measure=measure, stat=args.stat, group=label,
                                     orf=o, n_peptides=int(m.sum()),
                                     peak=float(np.nanmax(v[m])),
                                     median=float(np.nanmedian(v[m]))))

    table = pd.DataFrame(rows)
    csv = out_dir / f"{prefix}_genome_summary.csv"
    table.to_csv(csv, index=False)
    print(f"Wrote {csv}")

    for measure in args.measure:
        block = table[(table["measure"] == measure) & (table["view"] == "all")]
        if block.empty:
            continue
        print(f"\nPeak {measure} per ORF — whole library:")
        print(block.pivot_table(index="orf", columns="group",
                                values="peak").round(2).to_string())


if __name__ == "__main__":
    main()
