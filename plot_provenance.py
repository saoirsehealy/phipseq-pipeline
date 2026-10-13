#!/usr/bin/env python3
"""One-line provenance stamp for figures.

Several analyses in this project look similar but rest on very different inputs:
peptides chosen by Emma's ranking pipeline versus every peptide in the library,
and z-scores computed here versus phipflow's own data_zscore.csv. Those choices
change the numbers substantially, so every figure states which it used.

The stamp is width-aware. A fixed font size ran off the edge of narrow figures,
so each line is wrapped to the actual figure width and the font shrinks (to a
floor) before wrapping kicks in. `stamp` returns the y-coordinate of the top of
the block, so callers can place other footer text above it without colliding.
"""
import textwrap

FONT_MAX, FONT_MIN = 7.0, 5.4
CHAR_W = 0.60          # monospace advance width, in em

SRC = {
    "raw_cpm": "raw phipflow CPM (exp23/combined/wide_data/data_cpm_final_annotated.csv) | ALL matching peptides, no panel selection",
    "panel":   "Emma pipeline v1.7.2 unsupervised run | peptides = ranked/tuned panel from <prefix>_selected_panel.txt",
}
TRF = {
    "my_z":       "Z: computed in this repo, per-condition mocks, leave-one-out for mocks | NOT phipflow data_zscore.csv",
    "phipflow_z": "Z: phipflow data_zscore.csv (phippery binned beads-only model)",
    "log_cpm":    "transform: log10(CPM+1), no background correction",
    "cpm_minus":  "transform: CPM minus that peptide's max mock well (raw CPM scale)",
    "none":       "transform: none, raw CPM throughout",
    "rank_auc":   "AUC is rank-based, so identical for raw CPM and any per-peptide rescaling (asserted at runtime)",
}


def stamp(fig, source, transform, extra=None, y=0.004, x=0.004):
    """Add a provenance footer; return the y just above the block."""
    raw = [("SOURCE", SRC.get(source, source)),
           ("METHOD", TRF.get(transform, transform))]
    if extra:
        raw.append(("NOTE", extra))

    fig_w_pt = fig.get_figwidth() * 72 * (1 - 2 * x)
    fs = FONT_MAX
    while fs > FONT_MIN:
        if max(len(f"{k}  {v}") for k, v in raw) * fs * CHAR_W <= fig_w_pt:
            break
        fs -= 0.2
    width = max(24, int(fig_w_pt / (fs * CHAR_W)))

    lines = []
    for k, v in raw:
        wrapped = textwrap.wrap(f"{k}  {v}", width,
                                subsequent_indent=" " * (len(k) + 2))
        lines.extend(wrapped or [f"{k}  {v}"])
    txt = "\n".join(lines)
    fig.text(x, y, txt, ha="left", va="bottom", family="monospace",
             fontsize=fs, color="#6b6864", linespacing=1.45)

    # height of the block in figure fraction, plus a small gap
    line_h = fs * 1.45 / (fig.get_figheight() * 72)
    return y + len(lines) * line_h + 0.012
