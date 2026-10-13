#!/usr/bin/env python3
"""
PhIP-seq seroprevalence visualization

Generates four plot types per virus:
  1. Coverage plot       — max z-score per group as horizontal bar segments
  2. Seroprevalence      — % samples with z >= threshold per peptide, faceted by group
  3. Strain-exclusive    — seroprevalence for strain-exclusive peptides, faceted by strain
  4. Max signal per sample — boxplot of max z-score per sample, by group

Usage:
    python3 plot_phipseq_seroprevalence.py --exp-dir exp23 --viruses cchfv,ebola,hcv,rvfv

Input:  expXX/{virus}/{virus}[_{seg}]_relabelled.csv
Output: expXX/plots/{virus}_coverage.png
                    {virus}_seroprevalence.png
                    {virus}_strain_exclusive.png
                    {virus}_max_signal.png
"""

import argparse, re, zlib
from pathlib import Path
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
import seaborn as sns

# =============================================================================
# SAMPLE GROUP CONFIGURATION
# (regex_pattern, display_label, hex_color, include_in_plots)
# First match wins. Mock/Input excluded from seroprevalence plots.
# =============================================================================
SAMPLE_GROUPS = [
    (r'mock',    'Mock',      '#888888', False),
    (r'input',   'Input',     '#444444', False),
    (r'_neg\d',  'Negative',  '#AAAAAA', True),
    (r'CCHFV',   'CCHFV+',   '#E55934', True),
    (r'EBOV',    'EBOV+',    '#E8871A', True),
    (r'HCV',     'HCV+',     '#3A86FF', True),
    (r'RVFV',    'RVFV+',   '#8B5CF6', True),
]

ZSCORE_THRESHOLD = 3.5   # significance threshold for seroprevalence
ZSCORE_CAP       = None  # cap for coverage plot y-axis; None = no cap (autoscale)
MIN_GROUP_N      = 3     # minimum samples per group to plot seroprevalence

CPM_THRESHOLD    = 1.0   # log10(CPM+1) significance threshold for CPM-based plots (~CPM 9)
CPM_LABEL        = 'log10(CPM+1)'

# Auto-select focus group for strain-exclusive plot (keyed on virus name)
VIRUS_FOCUS_GROUP = {
    'ebola': 'EBOV+',
    'cchfv': 'CCHFV+',
    'hcv':   'HCV+',
    'rvfv':  'RVFV+',
}

# Default groups shown per virus (positive group + Negative control).
# Pass --all-groups to override and show every group.
VIRUS_SHOW_GROUPS = {
    'ebola': ['EBOV+', 'Negative'],
    'cchfv': ['CCHFV+', 'Negative'],
    'hcv':   ['HCV+', 'Negative'],
    'rvfv':  ['RVFV+', 'Negative'],
}


# =============================================================================
# DETERMINISTIC DOT JITTER
# =============================================================================

JITTER_SEED = 0   # change only if you want every dot layout to shift

def jitter(n, width, key=''):
    """Reproducible horizontal jitter for strip/box overlay dots.

    These plots previously drew jitter from the global np.random stream, which
    is never seeded — so dot positions depended on how many random draws had
    already happened in the process. Re-running a script, or generating one
    plot on its own instead of as part of a batch, moved every dot.

    Here the offsets come from a local RNG seeded by JITTER_SEED plus a stable
    checksum of `key`. The result depends only on (key, n), so a given group's
    dots land identically no matter which script drew them or what ran first.
    `key` should identify the plot and group (e.g. 'max_z|CCHFV+') so different
    panels don't all reuse one offset pattern.

    zlib.crc32 is used rather than hash() because hash() on str is salted per
    process (PYTHONHASHSEED), which would defeat the whole point.
    """
    rng = np.random.default_rng(JITTER_SEED + zlib.crc32(str(key).encode()))
    return rng.uniform(-width, width, size=n)


# =============================================================================
# SAMPLE GROUP HELPERS
# =============================================================================

def get_group(sample_name):
    for pattern, label, color, include in SAMPLE_GROUPS:
        if re.search(pattern, sample_name, re.IGNORECASE):
            return label, color, include
    return 'Other', '#CCCCCC', False   # exclude unmatched samples (controls etc.)


def build_group_info(columns):
    return {col: get_group(col) for col in columns}


# =============================================================================
# PEPTIDE NAME PARSING
# =============================================================================

def parse_relabelled_name(name):
    """
    Parse a genotyping-relabelled peptide name.
    Format (after original protein name):
        |lib_start-lib_end|virus|segment|strain_field[|strain2...][|blastp_coord]

    Returns dict with keys:
        lib_start, lib_end     : library peptide tile positions (aa)
        coord_start, coord_end : blastp-mapped positions in reference protein
        strain                 : strain label (or 'NoGenotype' / 'Ambiguous')
        exclusive              : True if mapped to exactly one strain
    """
    lib_m = re.search(r'\|(\d+)-(\d+)\|', name)
    if not lib_m:
        return None

    lib_start = int(lib_m.group(1))
    lib_end   = int(lib_m.group(2))
    remainder = name[lib_m.end():]
    parts     = remainder.split('|')   # [virus, segment, strain..., blastp_coord]

    if len(parts) < 3:
        return None

    blastp_raw = parts[-1]
    bm = re.match(r'^(\d+)-(\d+)$', blastp_raw)
    coord_start = int(bm.group(1)) if bm else None
    coord_end   = int(bm.group(2)) if bm else None

    strain_parts = parts[2:-1]   # everything between segment and blastp coord

    if not strain_parts or strain_parts[0] in ('NoGenotype', 'NOCOORD'):
        strain    = 'NoGenotype'
        exclusive = False
    elif strain_parts[0].startswith('Ambiguous'):
        strain    = 'Ambiguous'
        exclusive = False
    else:
        strain    = strain_parts[0]
        exclusive = True

    return {
        'lib_start':   lib_start,
        'lib_end':     lib_end,
        'coord_start': coord_start,
        'coord_end':   coord_end,
        'strain':      strain,
        'exclusive':   exclusive,
    }


def simplify_strain(strain_id):
    """Convert long accession-based strain ID to short readable label."""
    if not strain_id or strain_id in ('NoGenotype', 'Ambiguous', 'NOCOORD'):
        return strain_id
    # Remove leading accession prefix like NC_002549.1_ or DQ211646_
    s = re.sub(r'^[A-Z]{1,4}_?\d[\d\.]*_', '', strain_id)
    # Keep underscores for readability
    return s.strip('_')


# =============================================================================
# Z-SCORE COMPUTATION
# =============================================================================

def compute_zscore(df, mock_regex=r'mock'):
    """Per-peptide z-score: (cpm - mean_mock) / std_mock."""
    mock_cols = [c for c in df.columns if re.search(mock_regex, c, re.IGNORECASE)]
    if not mock_cols:
        print("  WARNING: no mock columns found — using all-sample mean")
        mock_cols = list(df.columns)
    mock_mean = df[mock_cols].mean(axis=1)
    mock_std  = df[mock_cols].std(axis=1)
    floor = mock_std[mock_std > 1e-6].median()
    mock_std[mock_std < 1e-6] = max(floor, 1e-6)
    return df.subtract(mock_mean, axis=0).divide(mock_std, axis=0)


def compute_log10_cpm(df):
    """Per-peptide log10(CPM+1). No mock background subtraction."""
    return np.log10(df.clip(lower=0) + 1)


# =============================================================================
# SEROPREVALENCE + CONFIDENCE INTERVALS
# =============================================================================

def wilson_ci(k, n, z_crit=1.96):
    """Wilson score interval. Returns (proportion_pct, lo_pct, hi_pct)."""
    if n == 0:
        return 0.0, 0.0, 0.0
    p = k / n
    denom  = 1 + z_crit**2 / n
    center = (p + z_crit**2 / (2 * n)) / denom
    margin = z_crit * np.sqrt(p * (1 - p) / n + z_crit**2 / (4 * n**2)) / denom
    return p * 100, max(0, center - margin) * 100, min(100, center + margin) * 100


def seroprevalence_per_peptide(z_group, threshold=ZSCORE_THRESHOLD):
    """
    For each peptide (row), compute % of samples (cols) with z >= threshold.
    Returns DataFrame with columns: sero, ci_lo, ci_hi, n
    """
    n = z_group.shape[1]
    rows = []
    for pep, vals in z_group.iterrows():
        k = int((vals >= threshold).sum())
        p, lo, hi = wilson_ci(k, n)
        rows.append({'peptide': pep, 'sero': p, 'ci_lo': lo, 'ci_hi': hi, 'n': n, 'k': k})
    return pd.DataFrame(rows).set_index('peptide')


# =============================================================================
# PLOT 1 — COVERAGE PLOT
# =============================================================================

def plot_coverage(z, group_info, meta, virus, out_path, stat='max',
                  metric_label='Z-Score', background_note=' (mock background)',
                  y_cap=ZSCORE_CAP):
    """
    Horizontal bar segments per peptide, colored by group.
    stat='max'    → bar height = max value across samples in that group
    stat='median' → bar height = median value across samples in that group
    """
    active = [(g, c) for g, c, inc in
              [group_info[col] for col in z.columns] if inc]
    group_color = {}
    for g, c in active:
        group_color.setdefault(g, c)

    fig, ax = plt.subplots(figsize=(14, 5))

    for group, color in group_color.items():
        cols = [c for c in z.columns if group_info[c][0] == group]
        if not cols:
            continue
        group_z = z[cols]

        for pep in group_z.index:
            info = meta.get(pep)
            if not info:
                continue
            x0 = info['coord_start'] if info['coord_start'] is not None else info['lib_start']
            x1 = info['coord_end']   if info['coord_end']   is not None else info['lib_end']
            vals = group_z.loc[pep]
            y = float(vals.max() if stat == 'max' else vals.median())
            if y_cap is not None:
                y = min(y, y_cap)
            if y <= 0:
                continue
            ax.plot([x0, x1], [y, y], color=color, lw=2.5, alpha=0.55, solid_capstyle='butt')

    patches = [mpatches.Patch(color=c, label=g) for g, c in group_color.items()]
    ax.legend(handles=patches, loc='upper right', fontsize=9, frameon=False, title='Sample Type')

    stat_label = 'Max' if stat == 'max' else 'Median'
    ax.set_xlabel("Amino Acid Position", fontsize=11)
    ax.set_ylabel(f"{stat_label} {metric_label} across samples{background_note}", fontsize=11)
    ax.set_title(f"{virus}: Antibody Enrichment Across Protein ({stat_label})", fontsize=13)
    ax.set_ylim(bottom=0)
    sns.despine(ax=ax)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  → {out_path.name}")


# =============================================================================
# PLOT 2 — SEROPREVALENCE BAR CHART (per group)
# =============================================================================

def plot_seroprevalence(z, group_info, meta, virus, out_path,
                        threshold=ZSCORE_THRESHOLD, metric_label='Z-score',
                        background_note=' | Mock SD background'):
    """
    Proportion of individuals with value >= threshold per peptide, faceted by group.
    """
    active_groups = [
        (g, c) for g, c, inc in
        [(g, c, inc) for g, c, inc in
         [group_info[col] for col in z.columns] if inc]
        if True
    ]
    # Deduplicate
    seen = {}
    for g, c in active_groups:
        seen.setdefault(g, c)

    groups_to_plot = [(g, c) for g, c in seen.items()
                      if sum(1 for col in z.columns if group_info[col][0] == g) >= MIN_GROUP_N]
    if not groups_to_plot:
        print("  No groups with enough samples to plot seroprevalence")
        return

    n_groups = len(groups_to_plot)
    fig, axes = plt.subplots(n_groups, 1, figsize=(14, 3.5 * n_groups),
                             sharex=False, sharey=True)
    if n_groups == 1:
        axes = [axes]

    # Get coordinates and spans for x-axis
    coords = {}
    spans  = {}
    for pep in z.index:
        info = meta.get(pep)
        if info:
            x0 = info['coord_start'] if info['coord_start'] is not None else info['lib_start']
            x1 = info['coord_end']   if info['coord_end']   is not None else info['lib_end']
            coords[pep] = x0
            spans[pep]  = max(5, x1 - x0)   # minimum 5 aa so bars are always visible

    for ax, (group, color) in zip(axes, groups_to_plot):
        cols = [c for c in z.columns if group_info[c][0] == group]
        n = len(cols)
        sero_df = seroprevalence_per_peptide(z[cols], threshold)

        xs, ws, ys, yerr_lo, yerr_hi = [], [], [], [], []
        for pep, row in sero_df.iterrows():
            x = coords.get(pep)
            if x is None:
                continue
            xs.append(x)
            ws.append(spans.get(pep, 25))
            ys.append(row['sero'])
            yerr_lo.append(row['sero'] - row['ci_lo'])
            yerr_hi.append(row['ci_hi'] - row['sero'])

        if xs:
            order = np.argsort(xs)
            xs      = np.array(xs)[order]
            ws      = np.array(ws)[order]
            ys      = np.array(ys)[order]
            yerr_lo = np.array(yerr_lo)[order]
            yerr_hi = np.array(yerr_hi)[order]

            ax.bar(xs, ys, width=ws, color=color, alpha=0.85, align='edge')
            ax.errorbar(xs + ws / 2, ys,
                        yerr=[np.clip(yerr_lo, 0, None), np.clip(yerr_hi, 0, None)],
                        fmt='none', color='black', lw=0.8, capsize=2)

        ax.set_ylim(0, 105)
        ax.yaxis.set_major_formatter(mticker.PercentFormatter())
        ax.set_ylabel("% with significant\nresponse", fontsize=9)
        ax.set_title(group, fontsize=10)
        ax.axhline(0, color='black', lw=0.5)
        sns.despine(ax=ax)

    axes[-1].set_xlabel("Amino Acid Position (reference protein)", fontsize=11)
    fig.suptitle(
        f"{virus}: Proportion of individuals with significant response\n"
        f"{metric_label} ≥ {threshold} threshold{background_note}",
        fontsize=12, y=1.01
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  → {out_path.name}")


# =============================================================================
# PLOT 3 — STRAIN-EXCLUSIVE SEROPREVALENCE
# =============================================================================

def plot_strain_exclusive(z, group_info, meta, virus, out_path,
                          threshold=ZSCORE_THRESHOLD, focus_group=None,
                          metric_label='Z-score'):
    """
    Seroprevalence for strain-exclusive peptides only, faceted by strain.
    focus_group: if set, only use samples from this group (e.g. 'EBOV+').
                 If None, uses the largest seropositive group.
    """
    # Pick which sample group to use (first with enough samples)
    active_groups = []
    for g, c, inc in [group_info[col] for col in z.columns]:
        if inc and g not in ('Negative',) and g not in [x[0] for x in active_groups]:
            n = sum(1 for col in z.columns if group_info[col][0] == g)
            if n >= MIN_GROUP_N:
                active_groups.append((g, c, n))

    if focus_group:
        chosen = [(g, c, n) for g, c, n in active_groups if g == focus_group]
    else:
        # Use the non-negative group with the most samples
        chosen = [max(
            [(g, c, n) for g, c, n in active_groups if g != 'Negative'],
            key=lambda x: x[2]
        )] if active_groups else []

    if not chosen:
        print("  No suitable group for strain-exclusive plot")
        return

    use_group, use_color, _ = chosen[0]
    use_cols = [c for c in z.columns if group_info[c][0] == use_group]
    z_group  = z[use_cols]

    # Collect exclusive peptides and group by strain
    strain_peptides = defaultdict(list)
    for pep in z_group.index:
        info = meta.get(pep)
        if info and info['exclusive'] and info['strain'] not in ('NoGenotype', None):
            strain = simplify_strain(info['strain'])
            strain_peptides[strain].append(pep)

    if not strain_peptides:
        print("  No strain-exclusive peptides found")
        return

    # Sort strains by number of exclusive peptides (most first)
    strains = sorted(strain_peptides.keys(), key=lambda s: -len(strain_peptides[s]))
    n_strains = len(strains)
    counts = {s: len(strain_peptides[s]) for s in strains}

    fig, axes = plt.subplots(n_strains, 1, figsize=(14, 3.2 * n_strains),
                             sharex=False, sharey=True)
    if n_strains == 1:
        axes = [axes]

    for ax, strain in zip(axes, strains):
        peps = strain_peptides[strain]
        sero_df = seroprevalence_per_peptide(z_group.loc[peps], threshold)

        xs, ws, ys, yerr_lo, yerr_hi = [], [], [], [], []
        for pep, row in sero_df.iterrows():
            info = meta[pep]
            x0 = info['coord_start'] if info['coord_start'] is not None else info['lib_start']
            x1 = info['coord_end']   if info['coord_end']   is not None else info['lib_end']
            xs.append(x0)
            ws.append(max(5, x1 - x0))
            ys.append(row['sero'])
            yerr_lo.append(row['sero'] - row['ci_lo'])
            yerr_hi.append(row['ci_hi'] - row['sero'])

        if xs:
            order = np.argsort(xs)
            xs      = np.array(xs)[order]
            ws      = np.array(ws)[order]
            ys      = np.array(ys)[order]
            yerr_lo = np.array(yerr_lo)[order]
            yerr_hi = np.array(yerr_hi)[order]

            ax.bar(xs, ys, width=ws, color=use_color, alpha=0.85, align='edge')
            ax.errorbar(xs + ws / 2, ys,
                        yerr=[np.clip(yerr_lo, 0, None), np.clip(yerr_hi, 0, None)],
                        fmt='none', color='black', lw=0.8, capsize=2)

        ax.set_ylim(0, 105)
        ax.yaxis.set_major_formatter(mticker.PercentFormatter())
        ax.set_ylabel("% AV samples\nseropositive", fontsize=9)
        ax.set_title(f"{strain}  (n={counts[strain]} exclusive peptides)", fontsize=10)
        ax.axhline(0, color='black', lw=0.5)
        sns.despine(ax=ax)

    axes[-1].set_xlabel(f"Amino acid position (strain-specific coordinates)", fontsize=11)
    strain_summary = ", ".join(f"{s}={counts[s]}" for s in strains[:8])
    fig.suptitle(
        f"{use_group} seroprevalence — strain-exclusive {virus} peptides\n"
        f"{metric_label} ≥ {threshold} | Exclusive peptides: {strain_summary}",
        fontsize=11, y=1.01
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  → {out_path.name}")


# =============================================================================
# PLOT 4 — PER-SAMPLE SIGNAL SUMMARY (median z + breadth)
# =============================================================================

def plot_signal(z, group_info, virus, out_path,
                threshold=ZSCORE_THRESHOLD, metric_label='Z'):
    """
    Two-panel violin + strip plot per sample group:
      Left  — median value across all virus peptides per sample
      Right — number of peptides with value >= threshold per sample (breadth)
    """
    records = []
    for col in z.columns:
        grp, color, include = group_info[col]
        if not include:
            continue
        vals = z[col]
        records.append({
            'sample':   col,
            'group':    grp,
            'color':    color,
            'median_z': float(vals.median()),
            'breadth':  int((vals >= threshold).sum()),
        })
    if not records:
        return

    df_plot = pd.DataFrame(records)
    actual_groups = set(df_plot['group'])
    group_order = [lbl for _, lbl, _, inc in SAMPLE_GROUPS if inc and lbl in actual_groups]
    if 'Mock' in actual_groups and 'Mock' not in group_order:
        group_order.append('Mock')
    if not group_order:
        group_order = sorted(actual_groups)
    palette = {r['group']: r['color']
               for _, r in df_plot[['group', 'color']].drop_duplicates().iterrows()}

    fig, ax = plt.subplots(figsize=(max(5, len(group_order) * 1.8), 5))

    box_data = [df_plot[df_plot['group'] == g]['breadth'].dropna().values
                for g in group_order]
    bp = ax.boxplot(box_data, positions=range(len(group_order)),
                    widths=0.4, patch_artist=True,
                    flierprops=dict(marker=''), medianprops=dict(lw=2, color='black'))
    for patch, g in zip(bp['boxes'], group_order):
        patch.set_facecolor(palette[g])
        patch.set_alpha(0.5)
    for g, x_pos in zip(group_order, range(len(group_order))):
        vals = df_plot[df_plot['group'] == g]['breadth'].values
        jx = x_pos + jitter(len(vals), 0.18, key=f'breadth|{g}')
        ax.scatter(jx, vals, color=palette[g], s=28, alpha=0.85, zorder=3)
    ax.set_xticks(range(len(group_order)))
    ax.set_xticklabels(group_order, rotation=20, ha='right', fontsize=9)
    ax.set_ylabel(f"Peptides with {metric_label} ≥ {threshold}", fontsize=11)
    ax.set_xlabel("")
    sns.despine(ax=ax)

    ax.set_title(f"{virus}: Breadth of response per sample", fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  → {out_path.name}")


# =============================================================================
# PLOT 5 — MAX SIGNAL PER SAMPLE
# =============================================================================

def plot_max_signal(z, group_info, virus, out_path,
                    threshold=ZSCORE_THRESHOLD, metric_label='Z-Score',
                    background_note=' (mock background)', stat='max'):
    """Per-sample summary across peptides, grouped by sample type — box + strip plot.

    `stat` selects the summary: 'max' (default, unchanged behaviour) reports each
    sample's single highest-scoring peptide; 'median' reports the typical peptide
    instead, which is far less sensitive to one noisy well.
    """
    if stat not in ('max', 'median'):
        raise ValueError(f"stat must be 'max' or 'median', got {stat!r}")
    summarise = (lambda s: s.max()) if stat == 'max' else (lambda s: s.median())
    records = []
    for col in z.columns:
        grp, color, include = group_info[col]
        if not include:
            continue
        records.append({
            'sample': col,
            'group':  grp,
            'color':  color,
            'max_z':  float(summarise(z[col])),
        })
    if not records:
        return
    df_plot = pd.DataFrame(records)

    actual_groups = set(df_plot['group'].values)
    group_order = [lbl for _, lbl, _, inc in SAMPLE_GROUPS
                   if inc and lbl in actual_groups]
    if 'Mock' in actual_groups and 'Mock' not in group_order:
        group_order.append('Mock')
    if not group_order:
        group_order = sorted(actual_groups)
    palette = {r['group']: r['color'] for _, r in df_plot[['group', 'color']].drop_duplicates().iterrows()}

    fig, ax = plt.subplots(figsize=(max(6, len(group_order) * 1.5), 5))
    # Use matplotlib boxplot to avoid seaborn 0.13 hue/boxprops bug
    box_data = [df_plot[df_plot['group'] == g]['max_z'].dropna().values
                for g in group_order]
    bp = ax.boxplot(box_data, positions=range(len(group_order)),
                    widths=0.45, patch_artist=True,
                    flierprops=dict(marker=''), medianprops=dict(lw=2))
    for patch, g in zip(bp['boxes'], group_order):
        patch.set_facecolor(palette[g])
        patch.set_alpha(0.5)
    for g, x_pos in zip(group_order, range(len(group_order))):
        vals = df_plot[df_plot['group'] == g]['max_z'].values
        jx = x_pos + jitter(len(vals), 0.18, key=f'{stat}_z|{g}' if stat != 'max' else f'max_z|{g}')
        ax.scatter(jx, vals, color=palette[g], s=30, alpha=0.85, zorder=3)
    ax.set_xticks(range(len(group_order)))
    ax.set_xticklabels(group_order)

    # Cap y-axis at 120% of the 99th percentile to prevent outliers squashing the main data
    all_vals = df_plot['max_z'].values
    y_cap = np.percentile(all_vals, 99) * 1.25
    y_cap = max(y_cap, threshold * 3)
    ax.set_ylim(bottom=min(0, np.percentile(all_vals, 1)), top=y_cap)

    ax.set_xlabel("")
    ax.set_ylabel(f"{stat.capitalize()} {metric_label}{background_note}", fontsize=11)
    ax.set_title(f"{virus}: {'Maximum' if stat == 'max' else 'Median'} antibody signal per sample", fontsize=13)
    ax.tick_params(axis='x', rotation=20)
    ax.axhline(threshold, color='grey', lw=1, ls='--', alpha=0.5,
               label=f'Threshold ({threshold})')
    ax.legend(fontsize=8, frameon=False)
    sns.despine(ax=ax)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  → {out_path.name}")


# =============================================================================
# TREATMENT / REPLICATE HELPERS
# =============================================================================

def parse_sample_trt(sample_name):
    """
    Return (trt_label, base_id, is_dup) from a sample name.
      trt1_EBOV1_AV571     → ('trt1', 'EBOV1_AV571', False)
      dup_trt2_neg1_HD9    → ('trt2', 'neg1_HD9',    True)
      trt3_mock1           → ('trt3', 'mock1',        False)
      inputOCCUSCAN_1      → (None,   name,            False)
    """
    name = sample_name
    is_dup = name.startswith('dup_')
    if is_dup:
        name = name[4:]
    m = re.match(r'trt(\d+)_(.*)', name, re.IGNORECASE)
    if not m:
        return None, sample_name, is_dup
    return f"trt{m.group(1)}", m.group(2), is_dup


def _build_breadth_df(z, group_info, allowed_groups=None,
                      skip_dup=True, threshold=ZSCORE_THRESHOLD):
    """
    Per-sample DataFrame: group, color, trt, base, breadth, max_z.
    allowed_groups: explicit set of group labels to include; None = use include flag.
    """
    records = []
    for col in z.columns:
        grp, color, include = group_info[col]
        if allowed_groups is not None:
            if grp not in allowed_groups:
                continue
        elif not include:
            continue
        trt, base, is_dup = parse_sample_trt(col)
        if skip_dup and is_dup:
            continue
        records.append({
            'sample':  col,
            'group':   grp,
            'color':   color,
            'trt':     trt,
            'base':    base,
            'breadth': int((z[col] >= threshold).sum()),
            'max_z':   float(z[col].max()),
        })
    return pd.DataFrame(records) if records else pd.DataFrame()


# =============================================================================
# PLOT 6 — REPLICATE CONCORDANCE (scatter trt1 vs trt2 vs trt3)
# =============================================================================

def plot_replicates(z, group_info, virus, out_path, metric='breadth',
                    threshold=ZSCORE_THRESHOLD, value_label='Z'):
    """
    Scatter of per-sample metric between every treatment pair.
    metric: 'breadth' (peptides value>=threshold) or 'max_z'.
    Each dot = one matched individual (same base name across trt).
    Diagonal = perfect concordance.
    """
    df = _build_breadth_df(z, group_info, threshold=threshold)
    if df.empty:
        print("  No data for replicate concordance plot")
        return

    trts = sorted(df['trt'].dropna().unique())
    if len(trts) < 2:
        print(f"  Only {len(trts)} treatment(s) — skipping replicate concordance")
        return

    pairs = [(trts[i], trts[j])
             for i in range(len(trts)) for j in range(i + 1, len(trts))]
    actual_groups = df['group'].unique()
    group_order = [lbl for _, lbl, _, inc in SAMPLE_GROUPS if inc and lbl in actual_groups]
    palette = {row['group']: row['color']
               for _, row in df[['group', 'color']].drop_duplicates().iterrows()}

    col = metric  # 'breadth' or 'max_z'
    metric_label = f"Peptides {value_label}≥{threshold}" if metric == 'breadth' else f"Max {value_label}"

    fig, axes = plt.subplots(1, len(pairs), figsize=(5 * len(pairs), 4.5), squeeze=False)

    for ax, (ta, tb) in zip(axes[0], pairs):
        sub_a = (df[df['trt'] == ta][['base', 'group', col]]
                 .rename(columns={col: 'va'}))
        sub_b = (df[df['trt'] == tb][['base', 'group', col]]
                 .rename(columns={col: 'vb'}))
        merged = sub_a.merge(sub_b, on=['base', 'group'])

        if merged.empty:
            ax.set_title(f"{ta} vs {tb}\n(no matched pairs)")
            continue

        for grp in group_order:
            pts = merged[merged['group'] == grp]
            if pts.empty:
                continue
            ax.scatter(pts['va'], pts['vb'],
                       color=palette.get(grp, '#AAAAAA'), label=grp,
                       s=55, alpha=0.85, zorder=3)

        vmax = max(merged['va'].max(), merged['vb'].max(), 1) * 1.1
        ax.plot([0, vmax], [0, vmax], 'k--', lw=0.8, alpha=0.35)
        ax.set_xlim(-1, vmax)
        ax.set_ylim(-1, vmax)
        ax.set_xlabel(f"{metric_label} ({ta})", fontsize=10)
        ax.set_ylabel(f"{metric_label} ({tb})", fontsize=10)
        ax.set_title(f"{ta} vs {tb}  (n={len(merged)} pairs)", fontsize=11)
        ax.legend(fontsize=8, frameon=False)
        sns.despine(ax=ax)

    fig.suptitle(
        f"{virus}: Replicate concordance — {metric_label.lower()}",
        fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  → {out_path.name}")


# =============================================================================
# PLOT 7 — CONDITION COMPARISON (trt1 / trt2 / trt3 per sample group)
# =============================================================================

def plot_conditions(z, group_info, virus, out_path, show_groups=None,
                    metric='breadth', threshold=ZSCORE_THRESHOLD, value_label='Z'):
    """
    Box+dot per-sample metric per treatment, faceted by sample group.
    Lines connect the same individual across conditions.
    metric: 'breadth' (peptides value>=threshold) or 'max_z'.
    show_groups: which group labels to include (e.g. ['EBOV+', 'Negative', 'Mock']).
    """
    df = _build_breadth_df(z, group_info, allowed_groups=show_groups, threshold=threshold)
    if df.empty:
        print("  No data for condition comparison")
        return

    df = df[df['trt'].notna()].copy()
    actual_groups = set(df['group'])

    # Build group order using labels (not patterns) from SAMPLE_GROUPS.
    # Order: virus-positive included groups, then Negative, then Mock.
    inc_groups  = [lbl for _, lbl, _, inc in SAMPLE_GROUPS
                   if inc and lbl in actual_groups and lbl != 'Negative']
    neg_groups  = ['Negative'] if 'Negative' in actual_groups else []
    mock_groups = [lbl for _, lbl, _, inc in SAMPLE_GROUPS
                   if not inc and lbl in actual_groups and lbl != 'Input']
    group_order = inc_groups + neg_groups + mock_groups
    if not group_order:
        group_order = sorted(actual_groups)

    palette = {row['group']: row['color']
               for _, row in df[['group', 'color']].drop_duplicates().iterrows()}

    col = metric
    metric_label = f"Peptides {value_label}≥{threshold}" if metric == 'breadth' else f"Max {value_label}"

    fig, axes = plt.subplots(1, len(group_order),
                             figsize=(max(8, 2.8 * len(group_order)), 5),
                             sharey=True, squeeze=False)

    for ax, grp in zip(axes[0], group_order):
        sub = df[df['group'] == grp].copy()
        grp_trts = sorted(sub['trt'].unique())
        if not grp_trts:
            ax.set_title(grp)
            continue

        box_data = [sub[sub['trt'] == t][col].values for t in grp_trts]
        bp = ax.boxplot(box_data, positions=range(len(grp_trts)),
                        widths=0.4, patch_artist=True,
                        flierprops=dict(marker=''),
                        medianprops=dict(lw=2, color='black'))
        for patch in bp['boxes']:
            patch.set_facecolor(palette.get(grp, '#AAAAAA'))
            patch.set_alpha(0.45)

        # Track each dot's actual jittered position for line drawing
        dot_pos = {}   # (base, trt) → (jx, y)
        for xi, t in enumerate(grp_trts):
            sub_t = sub[sub['trt'] == t]
            vals  = sub_t[col].values
            bases = sub_t['base'].values
            jx    = xi + jitter(len(vals), 0.15, key=f'cond|{grp}|{t}')
            ax.scatter(jx, vals, color=palette.get(grp, '#AAAAAA'),
                       s=30, alpha=0.85, zorder=3)
            for b, x, y in zip(bases, jx, vals):
                dot_pos[(b, t)] = (x, y)

        # Connect individuals present in ALL treatments using actual dot positions
        pivot = sub.pivot_table(
            index='base', columns='trt', values=col, aggfunc='first')
        complete = pivot.dropna(subset=grp_trts, how='any')
        line_color = palette.get(grp, '#AAAAAA')
        for base, _ in complete.iterrows():
            coords = [dot_pos.get((base, t)) for t in grp_trts]
            if any(c is None for c in coords):
                continue
            ax.plot([c[0] for c in coords], [c[1] for c in coords],
                    color=line_color, lw=0.7, alpha=0.45, zorder=2)

        ax.set_xticks(range(len(grp_trts)))
        ax.set_xticklabels(grp_trts, fontsize=9, rotation=20, ha='right')
        n_ind = sub['base'].nunique()
        ax.set_title(f"{grp}\n(n={n_ind})", fontsize=10)
        sns.despine(ax=ax)

    axes[0][0].set_ylabel(metric_label, fontsize=11)
    fig.suptitle(f"{virus}: Condition comparison — trt1 / trt2 / trt3  ({metric_label})",
                 fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  → {out_path.name}")


# =============================================================================
# MAIN
# =============================================================================

def find_relabelled_csv(exp_dir, virus):
    vdir = exp_dir / virus
    if not vdir.is_dir():
        return None
    candidates = sorted(vdir.glob(f"{virus}_*_relabelled.csv")) + \
                 [vdir / f"{virus}_relabelled.csv"]
    for c in candidates:
        if c.exists():
            return c
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--exp-dir',  required=True)
    ap.add_argument('--viruses',  default='cchfv,ebola,hcv,rvfv')
    ap.add_argument('--out-dir',  default=None)
    ap.add_argument('--threshold', type=float, default=ZSCORE_THRESHOLD,
                    help=f"Z-score threshold for seroprevalence (default {ZSCORE_THRESHOLD})")
    ap.add_argument('--focus-group', default=None,
                    help="Sample group to use for strain-exclusive plot (e.g. 'EBOV+')")
    ap.add_argument('--all-groups', action='store_true',
                    help="Show all sample groups in every plot (default: only virus-specific group + Negative)")
    ap.add_argument('--max-z-per-sample', type=float, default=None,
                    help="Exclude individual samples whose max z-score across virus peptides exceeds this value")
    args = ap.parse_args()

    exp_dir = Path(args.exp_dir).resolve()
    out_dir = Path(args.out_dir).resolve() if args.out_dir else exp_dir / 'plots'
    out_dir.mkdir(parents=True, exist_ok=True)
    viruses = [v.strip() for v in args.viruses.split(',')]

    for virus in viruses:
        print(f"\n{'='*55}")
        print(f"  {virus.upper()}")
        print(f"{'='*55}")

        csv_path = find_relabelled_csv(exp_dir, virus)
        if not csv_path:
            print(f"  ERROR: no relabelled CSV in {exp_dir / virus}/")
            continue

        print(f"  Loading {csv_path.name} ...")
        df = pd.read_csv(csv_path, index_col=0)
        print(f"  {df.shape[0]:,} peptides × {df.shape[1]} samples")

        group_info = build_group_info(df.columns)
        counts = Counter(g for g, _, _ in group_info.values())
        print(f"  Groups: {dict(counts)}")

        # Filter to virus peptides
        mask = df.index.str.contains(virus, case=False, na=False)
        df_v = df[mask].copy()
        print(f"  {len(df_v)} peptides match '{virus}'")
        if len(df_v) == 0:
            continue

        # Parse peptide metadata (coordinates, strain, exclusive)
        print("  Parsing peptide coordinates and strain assignments...")
        meta = {}
        for pep in df_v.index:
            info = parse_relabelled_name(pep)
            if info:
                meta[pep] = info

        coord_count = sum(1 for v in meta.values() if v['coord_start'] is not None)
        excl_count  = sum(1 for v in meta.values() if v.get('exclusive'))
        print(f"  {coord_count} peptides with blastp coordinates, {excl_count} strain-exclusive")

        # Sort by coordinate
        def sort_key(pep):
            info = meta.get(pep)
            if info and info['coord_start'] is not None:
                return info['coord_start']
            return 999999
        df_v = df_v.loc[sorted(df_v.index, key=sort_key)]

        # Compute z-scores
        print("  Computing z-scores from mock background...")
        z = compute_zscore(df_v)

        # Plottable samples: respect include flag, then optionally restrict to
        # virus-specific group + Negative (unless --all-groups requested)
        allowed_groups = None
        if not args.all_groups:
            allowed_groups = set(VIRUS_SHOW_GROUPS.get(virus.lower(), []))

        plot_cols = [
            c for c in z.columns
            if group_info[c][2]                                          # include=True
            and (allowed_groups is None or group_info[c][0] in allowed_groups)
        ]
        if not plot_cols:
            print("  WARNING: no plottable columns after group filter — falling back to all groups")
            plot_cols = [c for c in z.columns if group_info[c][2]]

        z_plot  = z[plot_cols]
        gi_plot = {c: group_info[c] for c in plot_cols}

        # Optionally exclude high-z outliers from the Negative group only
        # (e.g. misclassified or cross-reactive negatives distorting the baseline)
        if args.max_z_per_sample is not None:
            before = len(plot_cols)
            plot_cols = [
                c for c in plot_cols
                if group_info[c][0] != 'Negative'           # always keep non-negative groups
                or float(z[c].max()) <= args.max_z_per_sample
            ]
            removed = before - len(plot_cols)
            if removed:
                print(f"  Excluded {removed} Negative sample(s) with max z > {args.max_z_per_sample}")

        shown = Counter(group_info[c][0] for c in plot_cols)
        print(f"  Plotting groups: {dict(shown)}")

        vname = virus.upper()
        thr   = args.threshold

        print("  [1/4] Coverage plots (max + median)...")
        plot_coverage(z_plot, gi_plot, meta, vname,
                      out_dir / f"{virus}_coverage_max.png",    stat='max')
        plot_coverage(z_plot, gi_plot, meta, vname,
                      out_dir / f"{virus}_coverage_median.png", stat='median')

        print("  [2/4] Seroprevalence bar chart...")
        plot_seroprevalence(z_plot, gi_plot, meta, vname,
                            out_dir / f"{virus}_seroprevalence.png", threshold=thr)

        print("  [3/4] Strain-exclusive seroprevalence...")
        focus_grp = args.focus_group or VIRUS_FOCUS_GROUP.get(virus.lower())
        plot_strain_exclusive(z_plot, gi_plot, meta, vname,
                              out_dir / f"{virus}_strain_exclusive.png",
                              threshold=thr, focus_group=focus_grp)

        print("  [4/5] Per-sample signal summary...")
        plot_signal(z_plot, gi_plot, vname,
                    out_dir / f"{virus}_signal_breadth.png")

        print("  [5/5] Max signal per sample...")
        plot_max_signal(z_plot, gi_plot, vname,
                        out_dir / f"{virus}_max_signal.png")

        # ── New comparison plots ────────────────────────────────────────────
        # Build mock-inclusive columns once (used by several plots below)
        mock_cols = [c for c in z.columns if group_info[c][0] == 'Mock']
        z_with_mock = z[plot_cols + mock_cols]
        gi_with_mock = dict(gi_plot)
        gi_with_mock.update(
            {c: ('Mock', '#888888', True) for c in mock_cols})
        cond_groups = set(VIRUS_SHOW_GROUPS.get(virus.lower(), [])) | {'Mock'}

        print("  [6/10] Replicate concordance — breadth...")
        plot_replicates(z_plot, gi_plot, vname,
                        out_dir / f"{virus}_replicates.png",
                        metric='breadth')

        print("  [7/10] Replicate concordance — max Z-score...")
        plot_replicates(z_plot, gi_plot, vname,
                        out_dir / f"{virus}_replicates_maxz.png",
                        metric='max_z')

        sample_groups = set(VIRUS_SHOW_GROUPS.get(virus.lower(), []))

        print("  [8/10] Condition comparison (samples) — breadth...")
        plot_conditions(z, group_info, vname,
                        out_dir / f"{virus}_conditions_breadth.png",
                        show_groups=sample_groups, metric='breadth')

        print("  [9/10] Condition comparison (samples) — max Z-score...")
        plot_conditions(z, group_info, vname,
                        out_dir / f"{virus}_conditions_max.png",
                        show_groups=sample_groups, metric='max_z')

        print("  Condition comparison (mock) — breadth & max Z-score...")
        plot_conditions(z, group_info, vname,
                        out_dir / f"{virus}_mock_conditions_breadth.png",
                        show_groups={'Mock'}, metric='breadth')
        plot_conditions(z, group_info, vname,
                        out_dir / f"{virus}_mock_conditions_z.png",
                        show_groups={'Mock'}, metric='max_z')

        print("  [10/10] Signal with mock — breadth & max Z-score...")
        plot_signal(z_with_mock, gi_with_mock, vname,
                    out_dir / f"{virus}_signal_mock_breadth.png")
        plot_max_signal(z_with_mock, gi_with_mock, vname,
                        out_dir / f"{virus}_max_signal_with_mock.png")

    print(f"\nAll plots saved to: {out_dir}/")


if __name__ == '__main__':
    main()
