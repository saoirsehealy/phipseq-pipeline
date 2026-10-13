#!/usr/bin/env python3
# Emma Thomson 19/03/26; updated to v1.6.4

"""
PhIP-seq merged analysis script v1.7.2
------------------------------------
Input:
  1) filtered virus-specific CPM matrix from Step 1
  2) optional virus/output prefix

This script:
  - ranks discriminating peptides
  - tunes peptide panel size and peptide CPM cut-off
  - calls sample seropositivity (or reactivity in unsupervised mode)
  - plots a heatmap using the selected peptide panel

Updates in v1.5.4:
  - mean_3_5sd terminology is now used instead of the ambiguous label 'strict'
  - mean_3_5sd, q100_plus_median, and hybrid sample-level calls are computed independently
  - ELISA summary now returns a table and includes Fisher exact tests for binary calls
  - duplicated sample-group parsing has been removed
  - unsupervised peptide ranking above mock thresholds is vectorised
  - multiple heatmap modes are available, including Q100 + mock median correction
  - default heatmap mode is now option 7
  - a positive-value CSV table is written for each heatmap mode, with peptides as rows and samples as columns

Outputs:
  <virus>_ranked_peptides.csv
  <virus>_selected_panel.txt
  <virus>_panel_tuning.csv
  <virus>_seropositivity_report.csv
  <virus>_heatmap_<mode>.png
  <virus>_heatmap_labels_<mode>.csv
  <virus>_heatmap_positive_values_<mode>.csv

New in v1.5.4:
  - supports --default and --strict non-interactive presets
  - accepts sample-name regexes on the command line (quote alternation patterns, e.g. 'epstein|EBV')
  - supports minimal or full output sets
  - presets now auto-accept recommended final defaults for panel size, CPM threshold, and sumCPM thresholds
  - optional unified positive-peptide threshold can be applied to both reporting calls and heatmap sample suppression
  - fixes missing call column creation before binary conversion
  - fixes undefined call_mode_choice in interactive heatmap prompting
"""

import argparse
import csv
import math
import os
import re
import sys
import warnings
from statistics import mean, median, stdev, StatisticsError

import numpy as np
import pandas as pd
from pandas.errors import PerformanceWarning

warnings.filterwarnings("ignore", category=PerformanceWarning)
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, mannwhitneyu, fisher_exact

NON_SAMPLE_TRAILING_COLS = {
    "max_mock_cpm",
    "mean_mock_cpm",
    "prevalence_pct_nonmock",
    "ambiguous_score",
    "source_peptide_count",
    "orig_peptides",
}


DERIVED_SAMPLE_SUFFIXES = (
    "__ambiguous_score",
)

DERIVED_SAMPLE_INFIXES = (
    "__mock__",
    "__max_mock_cpm",
    "__mean_mock_cpm",
    "__prevalence_pct_nonmock",
)

def is_true_sample_column(colname):
    col = str(colname)
    if col == "peptide" or col in NON_SAMPLE_TRAILING_COLS:
        return False
    if any(col.endswith(suf) for suf in DERIVED_SAMPLE_SUFFIXES):
        return False
    if any(infix in col for infix in DERIVED_SAMPLE_INFIXES):
        return False
    return True


STRICT_SD_MULTIPLIER = 3.5

DEFAULT_PRESET = {
    "common_virus_mode": False,
    "mock_regex": "mock",
    "pos_regex": None,
    "neg_regex": None,
    "candidate_sizes": [5, 10, 15, 20, 25, 30, 40, 50],
    "candidate_cpms": [5, 10, 15, 20, 30],
    "panel_choice": None,
    "cpm_choice": None,
    "call_mode": "both",
    "strict_hit_thr": 5,
    "strict_sum_thr": None,
    "q100med_hit_thr": 5,
    "q100med_sum_thr": None,
    "heatmap_mode": "mock_q100_plus_median",
    "include_mocks_in_heatmap": False,
    "min_positive_peptides_for_display": 5,
    "z_clip": STRICT_SD_MULTIPLIER,
}

STRICT_PRESET = {
    "common_virus_mode": False,
    "mock_regex": "mock",
    "pos_regex": None,
    "neg_regex": None,
    "candidate_sizes": [10, 15, 20, 25, 30, 40, 50],
    "candidate_cpms": [5, 10, 15, 20, 30],
    "panel_choice": None,
    "cpm_choice": None,
    "call_mode": "q100_plus_median",
    "strict_hit_thr": 6,
    "strict_sum_thr": None,
    "q100med_hit_thr": 6,
    "q100med_sum_thr": None,
    "heatmap_mode": "mock_q100_plus_median",
    "include_mocks_in_heatmap": False,
    "min_positive_peptides_for_display": 5,
    "z_clip": 4.0,
}


def get_preset(use_default=False, use_strict=False):
    if use_strict:
        return dict(STRICT_PRESET)
    if use_default:
        return dict(DEFAULT_PRESET)
    return None


def clean_vals(vals):
    return [v for v in vals if not pd.isna(v)]


def safe_mean(vals):
    vals = clean_vals(vals)
    return mean(vals) if vals else float("nan")


def safe_median(vals):
    vals = clean_vals(vals)
    return median(vals) if vals else float("nan")


def safe_sd(vals):
    vals = clean_vals(vals)
    if len(vals) < 2:
        return 0.0
    try:
        return stdev(vals)
    except StatisticsError:
        return 0.0

def safe_mad(vals):
    vals = clean_vals(vals)
    if not vals:
        return 0.0
    med = median(vals)
    return median([abs(v - med) for v in vals])


def robust_sigma_from_mad(vals):
    return 1.4826 * safe_mad(vals)


def safe_iqr(vals):
    vals = sorted(clean_vals(vals))
    if len(vals) < 2:
        return 0.0
    mid = len(vals) // 2
    if len(vals) % 2 == 0:
        lower = vals[:mid]
        upper = vals[mid:]
    else:
        lower = vals[:mid]
        upper = vals[mid + 1:]
    q1 = median(lower) if lower else vals[0]
    q3 = median(upper) if upper else vals[-1]
    return q3 - q1


def prompt_with_default(prompt, default):
    # Separate prompts cleanly so consecutive input() calls do not visually run together
    # on buffered terminals.
    print("", flush=True)
    sys.stdout.write(f"{prompt} [default: {default}]: ")
    sys.stdout.flush()
    value = input().strip()
    return value if value else str(default)


def yes_no_prompt(prompt, default="n"):
    suffix = "Y/n" if default.lower() == "y" else "y/N"
    while True:
        value = input(f"{prompt} [{suffix}]: ").strip().lower()
        if not value:
            return default.lower() == "y"
        if value in {"y", "yes"}:
            return True
        if value in {"n", "no"}:
            return False
        print("Please answer y or n.")


def prompt_or_preset(prompt, default, preset_value=None):
    if preset_value is not None:
        return str(preset_value)
    return prompt_with_default(prompt, default)


def yes_no_or_preset(prompt, default="n", preset_value=None):
    if preset_value is not None:
        return bool(preset_value)
    return yes_no_prompt(prompt, default=default)


def parse_index_selection(selection_text, allowed_samples):
    selection_text = str(selection_text).strip()
    if not selection_text:
        return []
    chosen = []
    tokens = [t.strip() for t in selection_text.split(",") if t.strip()]
    for token in tokens:
        if "-" in token:
            start_s, end_s = token.split("-", 1)
            start_i = int(start_s)
            end_i = int(end_s)
            if start_i > end_i:
                start_i, end_i = end_i, start_i
            for idx in range(start_i, end_i + 1):
                if idx < 1 or idx > len(allowed_samples):
                    raise ValueError(f"Sample index out of range: {idx}")
                chosen.append(allowed_samples[idx - 1])
        else:
            idx = int(token)
            if idx < 1 or idx > len(allowed_samples):
                raise ValueError(f"Sample index out of range: {idx}")
            chosen.append(allowed_samples[idx - 1])
    return list(dict.fromkeys(chosen))


def print_indexed_samples(samples, title):
    print(f"\n{title}:")
    for i, sample in enumerate(samples, start=1):
        print(f"  {i:>3}. {sample}")


def sanitize_name(s):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s.strip().lower())


def derive_prefix(matrix_path, user_prefix=None):
    if user_prefix:
        return sanitize_name(user_prefix)
    base = os.path.basename(matrix_path)
    base = re.sub(r"_filtered_matrix\.csv$", "", base)
    base = re.sub(r"_all_cpm_matrix_filtered\.csv$", "", base)
    base = re.sub(r"_cpm_matrix_filtered\.csv$", "", base)
    base = re.sub(r"\.csv$", "", base)
    base = re.sub(r"^taxid_", "", base)
    base = re.sub(r"_all$", "", base)
    return sanitize_name(base)


def parse_search_terms(search_term):
    text = str(search_term).strip()
    if not text:
        return []
    # Backward-compatible summary terms only. Actual filtering is performed
    # by build_filter_mask(), which supports simple OR groups.
    parts = re.split(r"(?:\s+|_)*(?:\bAND\b|&)(?:\s+|_)*|[;,]+", text, flags=re.IGNORECASE)
    parts = [p.strip(" _()") for p in parts if p.strip(" _()")]
    return parts or [text]


def _strip_outer_parentheses(text):
    text = str(text).strip()
    while text.startswith("(") and text.endswith(")"):
        inner = text[1:-1].strip()
        depth = 0
        balanced = True
        for ch in inner:
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
                if depth < 0:
                    balanced = False
                    break
        if balanced and depth == 0:
            text = inner
        else:
            break
    return text


def _split_top_level_bool(text, operator):
    """Split a simple Boolean expression at top-level AND/OR tokens.

    This deliberately supports the expressions used for peptide filtering, for example:
      cchfv AND (|s| OR |m|) AND Africa3
    It is not intended to be a full Boolean parser.
    """
    text = str(text)
    op = operator.upper()
    parts = []
    buf = []
    depth = 0
    i = 0
    while i < len(text):
        ch = text[i]
        if ch == "(":
            depth += 1
            buf.append(ch)
            i += 1
            continue
        if ch == ")":
            depth = max(0, depth - 1)
            buf.append(ch)
            i += 1
            continue
        if depth == 0 and text[i:i+len(op)].upper() == op:
            before_ok = (i == 0) or text[i-1].isspace() or text[i-1] in "()"
            after_i = i + len(op)
            after_ok = (after_i >= len(text)) or text[after_i].isspace() or text[after_i] in "()"
            if before_ok and after_ok:
                part = "".join(buf).strip()
                if part:
                    parts.append(part)
                buf = []
                i = after_i
                continue
        buf.append(ch)
        i += 1
    part = "".join(buf).strip()
    if part:
        parts.append(part)
    return parts


def build_filter_mask(peptide_series, expression):
    """Build a peptide filter mask from simple AND/OR substring expressions.

    Supported examples:
      Africa3
      cchfv AND |m| AND Africa3
      cchfv AND (|s| OR |m|) AND Africa3
    Tokens are literal substrings, matched case-insensitively.
    """
    expression = str(expression or "").strip()
    if not expression:
        return pd.Series(True, index=peptide_series.index), []

    and_terms = _split_top_level_bool(expression, "AND")
    mask = pd.Series(True, index=peptide_series.index)
    parsed_terms = []

    for term in and_terms:
        term = _strip_outer_parentheses(term)
        or_terms = [_strip_outer_parentheses(t).strip() for t in _split_top_level_bool(term, "OR")]
        or_terms = [t for t in or_terms if t]
        if not or_terms:
            continue
        term_mask = pd.Series(False, index=peptide_series.index)
        for ot in or_terms:
            # Remove optional quote marks but keep pipes literally.
            token = ot.strip().strip('\"').strip("'")
            if not token:
                continue
            term_mask |= peptide_series.str.contains(re.escape(token), case=False, na=False)
        mask &= term_mask
        parsed_terms.append(or_terms if len(or_terms) > 1 else or_terms[0])

    return mask, parsed_terms


def _normalise_coordinate_pair(first, second):
    first = int(first)
    second = int(second)
    if first > second:
        first, second = second, first
    return first, second


def _extract_coordinate_from_text(text):
    s = str(text).strip()
    if not s:
        return None
    m = re.fullmatch(r"(\d+)\s*-\s*(\d+)", s)
    if m:
        return _normalise_coordinate_pair(m.group(1), m.group(2))
    m = re.search(r"(?<!\d)(\d+)\s*-\s*(\d+)(?!\d)", s)
    if m:
        return _normalise_coordinate_pair(m.group(1), m.group(2))
    return None


def extract_peptide_location(name):
    """Return a location label, preferring a valid terminal coordinate field.

    Primary behaviour:
    - use the final pipe-delimited field when it contains coordinates

    Fallback behaviour:
    - scan earlier pipe-delimited fields from right to left and use the first
      coordinate-bearing field found
    - if no pipe field contains coordinates, fall back to the original final
      field (or full name if no pipes are present)
    """
    s = str(name).strip()
    parts = [p.strip() for p in s.split("|")] if "|" in s else [s]

    if parts:
        terminal = parts[-1].strip()
        if _extract_coordinate_from_text(terminal) is not None:
            return terminal

    for part in reversed(parts[:-1]):
        if _extract_coordinate_from_text(part) is not None:
            return part.strip()

    if parts:
        return parts[-1].strip()
    return s


def location_sort_key(location):
    """Return a numeric sort key for peptide locations.

    Primary behaviour:
    - parse an explicit coordinate pair from the provided location string

    Fallback behaviour:
    - recover trailing numeric tokens as before
    - if still unavailable, return a large sentinel so invalid rows sort last
    """
    s = str(location).strip()
    coord = _extract_coordinate_from_text(s)
    if coord is not None:
        first, second = coord
        return (first, second, s)

    nums = re.findall(r"\d+", s)
    if nums:
        if len(nums) >= 2:
            first, second = _normalise_coordinate_pair(nums[-2], nums[-1])
        else:
            first = int(nums[-1])
            second = first
        return (first, second, s)
    return (10**12, 10**12, s)


def extract_ambiguous_score(name):
    m = re.search(r"Ambiguous(\d+)", str(name), flags=re.IGNORECASE)
    return int(m.group(1)) if m else 0


def collapse_matrix_to_locations(df_in, sample_cols, mock_cols=None, selection_mode="raw_cpm", selection_floor=5.0):
    """Collapse peptides sharing the same terminal location.

    Behaviour:
    - peptides are collapsed ONLY when their extracted terminal location is identical
    - for each sample column, choose one winning peptide within that location group
    - choice can be made on raw CPM or a mock-aware score
    - all linked metadata, including mock CPMs for every mock sample, are taken from
      that SAME winning peptide row
    """
    work = df_in.copy()
    work["peptide_location"] = work["peptide"].astype(str).map(extract_peptide_location)
    work["ambiguous_score"] = work["peptide"].astype(str).map(extract_ambiguous_score)

    sample_cols_present = [c for c in sample_cols if c in work.columns]
    mock_cols_present = [c for c in (mock_cols or []) if c in work.columns]
    metadata_cols = [c for c in ["max_mock_cpm", "mean_mock_cpm", "prevalence_pct_nonmock"] if c in work.columns]

    def choose_best_index(group, sample):
        vals = pd.to_numeric(group[sample], errors="coerce")
        if not vals.notna().any():
            return None

        if selection_mode == "raw_cpm" or not mock_cols_present:
            return vals.idxmax()

        mock_block = group[mock_cols_present].apply(pd.to_numeric, errors="coerce")
        mock_mean = mock_block.mean(axis=1)
        mock_median = mock_block.median(axis=1)
        mock_q75 = mock_block.quantile(0.75, axis=1)
        mock_q100 = mock_block.max(axis=1)
        mock_sd = mock_block.std(axis=1).fillna(0.0)

        if selection_mode == "minus_mean_floor":
            score = vals - mock_mean - float(selection_floor)
        elif selection_mode == "minus_q75_floor":
            score = vals - mock_q75 - float(selection_floor)
        elif selection_mode == "minus_q100_floor":
            score = vals - mock_q100 - float(selection_floor)
        elif selection_mode == "minus_q100_plus_median_floor":
            score = vals - mock_q100 - mock_median - float(selection_floor)
        elif selection_mode == "zscore":
            score = pd.Series(np.where(mock_sd > 0, (vals - mock_mean) / mock_sd, -np.inf), index=group.index)
        else:
            score = vals

        score = pd.to_numeric(score, errors="coerce")
        if not score.notna().any():
            return vals.idxmax()
        # break ties on raw CPM so we keep the strongest actual signal if mock-aware scores tie
        score_df = pd.DataFrame({"score": score, "raw": vals}, index=group.index)
        score_df = score_df.sort_values(["score", "raw"], ascending=[False, False])
        return score_df.index[0]

    collapsed_rows = []
    for loc, group in work.groupby("peptide_location", sort=False):
        row = {
            "peptide": loc,
            "orig_peptides": " || ".join(map(str, group["peptide"].tolist())),
            "source_peptide_count": int(len(group)),
        }

        per_sample_amb_cols = []
        for sample in sample_cols_present:
            per_sample_amb_col = f"{sample}__ambiguous_score"
            per_sample_amb_cols.append(per_sample_amb_col)
            idx = choose_best_index(group, sample)
            if idx is not None:
                best = group.loc[idx]
                row[sample] = best[sample]
                row[per_sample_amb_col] = pd.to_numeric(best.get("ambiguous_score", np.nan), errors="coerce")
                for meta in metadata_cols:
                    row[f"{sample}__{meta}"] = best.get(meta, np.nan)
                for mock in mock_cols_present:
                    row[f"{sample}__mock__{mock}"] = best.get(mock, np.nan)
            else:
                row[sample] = np.nan
                row[per_sample_amb_col] = np.nan
                for meta in metadata_cols:
                    row[f"{sample}__{meta}"] = np.nan
                for mock in mock_cols_present:
                    row[f"{sample}__mock__{mock}"] = np.nan

        amb_vals = [row[c] for c in per_sample_amb_cols if c in row and not pd.isna(row[c])]
        row["ambiguous_score"] = float(np.nanmedian(amb_vals)) if amb_vals else np.nan

        for meta in metadata_cols:
            per_sample_meta_cols = [f"{s}__{meta}" for s in sample_cols_present if f"{s}__{meta}" in row]
            vals = [row[c] for c in per_sample_meta_cols if not pd.isna(row[c])]
            row[meta] = float(np.nanmax(vals)) if vals else np.nan

        collapsed_rows.append(row)

    collapsed = pd.DataFrame(collapsed_rows)
    collapsed["location_order"] = collapsed["peptide"].map(location_sort_key)
    collapsed = collapsed.sort_values("location_order").drop(columns=["location_order"]).reset_index(drop=True)
    return collapsed


def get_linked_mock_values_from_row(row, sample_name, mock_cols):
    """
    Return mock CPM values linked to the SAME peptide row selected for this sample.
    If sample-linked mock columns are absent, fall back to the raw mock columns.
    """
    linked_cols = [f"{sample_name}__mock__{m}" for m in mock_cols]
    if any(c in row.index for c in linked_cols):
        vals = [row.get(c, np.nan) for c in linked_cols]
    else:
        vals = [row.get(m, np.nan) for m in mock_cols]
    return np.array([float(v) if not pd.isna(v) else np.nan for v in vals], dtype=float)


def get_linked_mock_stats_from_row(row, sample_name, mock_cols):
    mock_vals = get_linked_mock_values_from_row(row, sample_name, mock_cols)
    clean = mock_vals[~np.isnan(mock_vals)]
    if clean.size == 0:
        return {
            "vals": clean,
            "mean": 0.0,
            "sd": 0.0,
            "median": 0.0,
            "q75": 0.0,
            "q100": 0.0,
            "positive_median": 0.0,
        }
    pos_vals = clean[clean > 0]
    pos_med = float(np.median(pos_vals)) if pos_vals.size > 0 else float(np.mean(clean))
    return {
        "vals": clean,
        "mean": float(np.mean(clean)),
        "sd": float(np.std(clean, ddof=1)) if clean.size > 1 else 0.0,
        "median": float(np.median(clean)),
        "q75": float(np.percentile(clean, 75)),
        "q100": float(np.max(clean)),
        "positive_median": pos_med,
    }

def shorten_peptide_name(name):
    parts = str(name).split("_")
    if len(parts) >= 4:
        return "_".join(parts[:4])
    return str(name)


def make_unique(names):
    seen = {}
    out = []
    for name in names:
        if name not in seen:
            seen[name] = 1
            out.append(name)
        else:
            seen[name] += 1
            out.append(f"{name}_{seen[name]}")
    return out


def mean_3_5sd_threshold(vals, floor=5.0, sd_multiplier=STRICT_SD_MULTIPLIER):
    vals = clean_vals(vals)
    if not vals:
        return floor
    return max(floor, safe_mean(vals) + sd_multiplier * safe_sd(vals))


# Backward-compatible alias retained for older internal references.
strict_threshold = mean_3_5sd_threshold


def q100_plus_median_threshold(vals, floor=5.0):
    vals = clean_vals(vals)
    if not vals:
        return floor
    return max(floor, max(vals) + safe_median(vals))


def build_peptide_thresholds(panel, peptide_presence_cpm, threshold_mode):
    peptide_thresholds = {}
    for pep in panel:
        row = matrix_lookup[pep]
        mvals = [row[c] for c in mock_cols if not pd.isna(row[c])]
        if threshold_mode in {"mean_3_5sd", "strict"}:
            threshold = mean_3_5sd_threshold(mvals, floor=peptide_presence_cpm)
        elif threshold_mode == "q100_plus_median":
            threshold = q100_plus_median_threshold(mvals, floor=peptide_presence_cpm)
        else:
            raise ValueError(f"Unknown threshold mode: {threshold_mode}")
        peptide_thresholds[pep] = threshold
    return peptide_thresholds



def compute_per_sample_metrics(panel, peptide_presence_cpm, threshold_mode):
    if not panel:
        return []

    if threshold_mode in {"mean_3_5sd", "strict"}:
        hit_col = "n_above_mock_mean_3_5sd_panel"
        sum_hits_col = "sum_cpm_panel_hits_only_mean_3_5sd"
        sum_log_hits_col = "sum_log10_cpm_panel_hits_only_mean_3_5sd"
    else:
        hit_col = "n_above_mock_q100_plus_median_panel"
        sum_hits_col = "sum_cpm_panel_hits_only_q100_plus_median"
        sum_log_hits_col = "sum_log10_cpm_panel_hits_only_q100_plus_median"

    per_sample = []
    for s in sample_cols:
        n_hit = 0
        sum_cpm_panel = 0.0
        sum_cpm_hits_only = 0.0
        sum_log10_cpm_hits_only = 0.0
        for pep in panel:
            row = matrix_lookup[pep]
            v = row.get(s, np.nan)
            if pd.isna(v):
                continue
            v = float(v)
            sum_cpm_panel += v

            mock_stats = get_linked_mock_stats_from_row(pd.Series(row), s, mock_cols)
            if threshold_mode in {"mean_3_5sd", "strict"}:
                thr = max(float(peptide_presence_cpm), mock_stats["mean"] + STRICT_SD_MULTIPLIER * mock_stats["sd"])
            else:
                thr = max(float(peptide_presence_cpm), mock_stats["q100"] + mock_stats["median"])

            if v > thr:
                n_hit += 1
                sum_cpm_hits_only += v
                sum_log10_cpm_hits_only += math.log10(v + 1.0)

        if s in mock_cols:
            group = "mock"
        elif s in pos_cols:
            group = "positive"
        elif s in neg_cols:
            group = "negative"
        else:
            group = "other"

        per_sample.append({
            "sample": s,
            "group": group,
            "n_panel_peptides": len(panel),
            hit_col: n_hit,
            "sum_cpm_panel": sum_cpm_panel,
            sum_hits_col: sum_cpm_hits_only,
            sum_log_hits_col: sum_log10_cpm_hits_only,
            "hit_fraction_panel": n_hit / len(panel) if panel else 0.0,
        })

    return per_sample


def derive_call_thresholds(per_sample, hit_col, sum_col, use_mocks_for_thresholds):
    if use_mocks_for_thresholds:
        threshold_source_label = "mocks"
        source_hits = [r[hit_col] for r in per_sample if r["group"] == "mock"]
        source_sum = [r[sum_col] for r in per_sample if r["group"] == "mock"]
    else:
        threshold_source_label = "negative comparison samples"
        source_hits = [r[hit_col] for r in per_sample if r["group"] == "negative"]
        source_sum = [r[sum_col] for r in per_sample if r["group"] == "negative"]

    src_hit_median = safe_median(source_hits)
    src_hit_sd = safe_sd(source_hits)
    src_hit_max = max(source_hits) if source_hits else 0

    src_sum_median = safe_median(source_sum)
    src_sum_sd = safe_sd(source_sum)
    src_sum_max = max(source_sum) if source_sum else 0.0

    hit_threshold_default = max(1, math.ceil(src_hit_median + STRICT_SD_MULTIPLIER * src_hit_sd))
    sum_threshold_default = src_sum_median + STRICT_SD_MULTIPLIER * src_sum_sd

    return {
        "threshold_source_label": threshold_source_label,
        "source_hits": source_hits,
        "source_sum": source_sum,
        "src_hit_median": src_hit_median,
        "src_hit_sd": src_hit_sd,
        "src_hit_max": src_hit_max,
        "src_sum_median": src_sum_median,
        "src_sum_sd": src_sum_sd,
        "src_sum_max": src_sum_max,
        "hit_threshold_default": hit_threshold_default,
        "sum_threshold_default": sum_threshold_default,
    }


def apply_calls_to_metrics(per_sample, hit_col, sum_col, hit_thr, sum_thr, unsupervised_mode):
    out = []
    for r in per_sample:
        rr = dict(r)
        hit_pass = rr[hit_col] >= hit_thr
        sum_pass = rr[sum_col] > sum_thr

        if unsupervised_mode:
            if hit_pass and sum_pass:
                call = "Reactive"
            elif rr[hit_col] == 0:
                call = "Non-reactive"
            else:
                call = "Indeterminate"
        else:
            if hit_pass and sum_pass:
                call = "Positive"
            elif rr[hit_col] == 0:
                call = "Negative"
            else:
                call = "Indeterminate"

        rr["call"] = call
        out.append(rr)
    return out


def binary_call_relaxed(call_value):
    call_value = str(call_value)
    if call_value in {"Positive", "Reactive", "Indeterminate"}:
        return 1
    return 0


def binary_call_strict(call_value):
    call_value = str(call_value)
    if call_value in {"Positive", "Reactive"}:
        return 1
    return 0


def derive_hybrid_call_thresholds(score_by_sample, hit_by_sample, neg_samples, all_nonmock_samples):
    neg_samples = [s for s in neg_samples if s in score_by_sample]
    if len(neg_samples) >= 1:
        source_label = "true negative controls"
        source_scores = [float(score_by_sample[s]) for s in neg_samples]
        med = safe_median(source_scores)
        sigma = robust_sigma_from_mad(source_scores)
        if sigma == 0.0:
            sigma = max(safe_sd(source_scores), 1e-9)
        t_low = med + 3.0 * sigma
        t_high = med + 5.0 * sigma
        return {
            "source_label": source_label,
            "source_scores": source_scores,
            "t_low": float(t_low),
            "t_high": float(t_high),
            "warning": " fewer than 5 true negative controls; thresholds may be unstable." if len(neg_samples) < 5 else ""
        }

    source_label = "all non-control samples (provisional fallback)"
    usable = [float(score_by_sample[s]) for s in all_nonmock_samples if s in score_by_sample]
    usable = sorted(clean_vals(usable))
    if not usable:
        return {
            "source_label": source_label,
            "source_scores": [],
            "t_low": 0.0,
            "t_high": 0.0,
            "warning": " no usable samples available for fallback thresholds."
        }
    q1 = float(np.percentile(usable, 25))
    q3 = float(np.percentile(usable, 75))
    iqr = q3 - q1
    t_low = q1 + 1.5 * iqr
    t_high = q1 + 2.5 * iqr
    return {
        "source_label": source_label,
        "source_scores": usable,
        "t_low": float(t_low),
        "t_high": float(t_high),
        "warning": " no true negative controls supplied; thresholds are provisional."
    }


def hybrid_call_from_score_and_hits(score, hit_count, t_low, t_high, unsupervised_mode=False,
                                    h_min_negative=3, h_min_positive=5):
    try:
        score = float(score)
    except Exception:
        score = 0.0
    try:
        hit_count = int(hit_count)
    except Exception:
        hit_count = 0

    if hit_count < h_min_negative or score < t_low:
        return "Non-reactive" if unsupervised_mode else "Negative"
    if hit_count >= h_min_positive and score >= t_high:
        return "Reactive" if unsupervised_mode else "Positive"
    return "Indeterminate"


def compute_mock_zscore_metrics(matrix_df, panel_peptides, sample_cols, mock_cols, cpm_floor):
    sub = matrix_df.set_index("peptide").loc[panel_peptides].copy()

    out = {
        "sum_cpm_panel_thresholded": {},
        "sum_z_mock_panel": {},
        "mean_z_mock_panel": {},
        "sum_z_mock_panel_thresholded": {},
        "mean_z_mock_panel_thresholded": {},
        "n_z_mock_above_2": {},
        "n_z_mock_above_3": {},
    }

    for s in sample_cols:
        sum_cpm_thr = 0.0
        z_vals_all = []
        z_vals_thr = []
        n_above_2 = 0
        n_above_3 = 0

        for pep in panel_peptides:
            row = sub.loc[pep]
            v = pd.to_numeric(row.get(s, np.nan), errors="coerce")
            if pd.isna(v):
                continue
            v = float(v)
            if v >= float(cpm_floor):
                sum_cpm_thr += v

            mock_stats = get_linked_mock_stats_from_row(row, s, mock_cols)
            mu = mock_stats["mean"]
            sd = mock_stats["sd"]
            if sd > 0:
                z = (v - mu) / sd
            else:
                z = 0.0
            if not np.isfinite(z):
                z = 0.0
            z_vals_all.append(float(z))

            if v >= float(cpm_floor):
                z_thr = float(z)
            else:
                z_thr = 0.0
            z_vals_thr.append(z_thr)
            if z_thr > 2:
                n_above_2 += 1
            if z_thr > 3:
                n_above_3 += 1

        out["sum_cpm_panel_thresholded"][s] = float(sum_cpm_thr)
        out["sum_z_mock_panel"][s] = float(np.sum(z_vals_all)) if z_vals_all else 0.0
        out["mean_z_mock_panel"][s] = float(np.mean(z_vals_all)) if z_vals_all else 0.0
        out["sum_z_mock_panel_thresholded"][s] = float(np.sum(z_vals_thr)) if z_vals_thr else 0.0
        out["mean_z_mock_panel_thresholded"][s] = float(np.mean(z_vals_thr)) if z_vals_thr else 0.0
        out["n_z_mock_above_2"][s] = int(n_above_2)
        out["n_z_mock_above_3"][s] = int(n_above_3)

    return out



def compute_hits_only_metrics(matrix_df, panel_peptides, sample_cols, mock_cols, cpm_floor, threshold_mode):
    sub = matrix_df.set_index("peptide").loc[panel_peptides].copy()

    if threshold_mode in {"mean_3_5sd", "strict"}:
        metric_suffix = "mean_3_5sd"
    else:
        metric_suffix = "q100_plus_median"

    sum_cpm_hits_only = {}
    sum_z_hits_only = {}
    n_hits_only = {}

    for s in sample_cols:
        total_cpm = 0.0
        total_z = 0.0
        n_hits = 0
        for pep in panel_peptides:
            row = sub.loc[pep]
            v = pd.to_numeric(row.get(s, np.nan), errors="coerce")
            if pd.isna(v):
                continue
            v = float(v)
            mock_stats = get_linked_mock_stats_from_row(row, s, mock_cols)
            if threshold_mode in {"mean_3_5sd", "strict"}:
                thr = max(float(cpm_floor), mock_stats["mean"] + 3.5 * mock_stats["sd"])
            else:
                thr = max(float(cpm_floor), mock_stats["q100"] + mock_stats["median"])

            sd = mock_stats["sd"]
            z = ((v - mock_stats["mean"]) / sd) if sd > 0 else 0.0
            if not np.isfinite(z):
                z = 0.0

            if v >= thr:
                n_hits += 1
                total_cpm += v
                total_z += z

        sum_cpm_hits_only[s] = float(total_cpm)
        sum_z_hits_only[s] = float(total_z)
        n_hits_only[s] = int(n_hits)

    return {
        f"sum_cpm_panel_hits_only_{metric_suffix}": sum_cpm_hits_only,
        f"sum_z_mock_hits_only_{metric_suffix}": sum_z_hits_only,
        f"n_hits_only_{metric_suffix}": n_hits_only,
    }



def per_row_positive_mock_median(mock_frame, threshold):
    baselines = []
    for _, row in mock_frame.iterrows():
        vals = row.astype(float).values
        pos_vals = vals[vals > threshold]
        if len(pos_vals) > 0:
            baselines.append(float(np.median(pos_vals)))
        else:
            baselines.append(float(np.mean(vals)) if len(vals) > 0 else 0.0)
    return np.array(baselines, dtype=float)



def compute_heatmap_transform_panel_metrics(matrix_df, panel_peptides, sample_cols, mock_cols, cpm_floor, z_clip=None):
    """Compute per-sample panel summaries for heatmap-style transforms.

    Key continuous corrected metric in v1.7.2:
      transform9_minus_mock_median_log10_sum =
      sum log10(max(0, CPM - per-peptide mock median) + 1)

    This uses the median mock CPM for the same peptide. When location collapse
    has been used, get_linked_mock_stats_from_row() supplies mock values from
    the same original peptide variant selected for that sample/location.
    """
    sub = matrix_df.set_index("peptide").loc[panel_peptides].copy()
    residual_floor = float(cpm_floor)
    if z_clip is None:
        z_clip = 4.0

    out = {
        "transform1_raw_log10_cpm_sum": {},
        "transform2_mock_zscore_sum": {},
        "transform3_minus_mock_mean_floor_log10_sum": {},
        "transform4_minus_positive_mock_median_floor_log10_sum": {},
        "transform5_minus_mock_q75_floor_log10_sum": {},
        "transform6_minus_mock_q100_floor_log10_sum": {},
        "transform7_minus_q100_plus_median_floor_log10_sum": {},
        "transform8_strict_mean3_5sd_only_log10_sum": {},
        "transform9_minus_mock_median_log10_sum": {},
    }

    for s in sample_cols:
        sums = {k: 0.0 for k in out.keys()}
        for pep in panel_peptides:
            row = sub.loc[pep]
            v = pd.to_numeric(row.get(s, np.nan), errors="coerce")
            if pd.isna(v):
                continue
            v = float(v)
            stats = get_linked_mock_stats_from_row(row, s, mock_cols)

            t1 = math.log10(v + 1.0)

            if stats["sd"] > 0:
                z = (v - stats["mean"]) / stats["sd"]
            else:
                z = z_clip if v > stats["mean"] else 0.0
            t2 = min(max(z, 0.0), z_clip)

            def baseline_val(baseline):
                return math.log10(max(0.0, v - baseline - residual_floor) + 1.0)

            t3 = baseline_val(stats["mean"])
            t4 = baseline_val(stats["positive_median"])
            t5 = baseline_val(stats["q75"])
            t6 = baseline_val(stats["q100"])
            t7 = baseline_val(stats["q100"] + stats["median"])
            strict_thr = max(float(cpm_floor), stats["mean"] + STRICT_SD_MULTIPLIER * stats["sd"])
            t8 = math.log10(v + 1.0) if v > strict_thr else 0.0
            # v1.7.2: biologically interpretable continuous corrected signal.
            # This subtracts the median mock CPM for this specific peptide only.
            # No q100/max penalty and no CPM-floor subtraction are applied here.
            t9 = math.log10(max(0.0, v - stats["median"]) + 1.0)

            sums["transform1_raw_log10_cpm_sum"] += t1
            sums["transform2_mock_zscore_sum"] += t2
            sums["transform3_minus_mock_mean_floor_log10_sum"] += t3
            sums["transform4_minus_positive_mock_median_floor_log10_sum"] += t4
            sums["transform5_minus_mock_q75_floor_log10_sum"] += t5
            sums["transform6_minus_mock_q100_floor_log10_sum"] += t6
            sums["transform7_minus_q100_plus_median_floor_log10_sum"] += t7
            sums["transform8_strict_mean3_5sd_only_log10_sum"] += t8
            sums["transform9_minus_mock_median_log10_sum"] += t9

        for k in out.keys():
            out[k][s] = float(sums[k])

    return out


def selected_continuous_signal_metric_name(method):
    """Return the continuous per-sample metric that corresponds to the selected
    heatmap/ranking method. For mock-corrected CPM methods this is a summed
    log10(max(0, CPM - mock-derived baseline - CPM floor) + 1) signal.
    """
    method = str(method or "mock_q100_plus_median").strip().lower()
    aliases = {
        "1": "raw_cpm",
        "2": "zscore",
        "3": "minus_mean",
        "4": "positive_mock_median",
        "5": "mock_q75",
        "6": "mock_q100",
        "7": "mock_q100_plus_median",
        "8": "mean_3_5sd",
        "9": "mock_median",
    }
    method = aliases.get(method, method)
    if method in {"minus_mean", "mock_mean"}:
        return "transform3_minus_mock_mean_floor_log10_sum"
    if method in {"mock_median", "median", "minus_median"}:
        return "transform9_minus_mock_median_log10_sum"
    if method in {"positive_mock_median", "mock_positive_median"}:
        return "transform4_minus_positive_mock_median_floor_log10_sum"
    if method == "mock_q75":
        return "transform5_minus_mock_q75_floor_log10_sum"
    if method == "mock_q100":
        return "transform6_minus_mock_q100_floor_log10_sum"
    if method == "mock_q100_plus_median":
        return "transform7_minus_q100_plus_median_floor_log10_sum"
    if method == "mean_3_5sd":
        return "transform8_strict_mean3_5sd_only_log10_sum"
    if method == "zscore":
        return "transform2_mock_zscore_sum"
    return "transform1_raw_log10_cpm_sum"

def normalise_sample_for_merge(sample_name):
    s = str(sample_name).strip()
    # remove common run prefixes such as exp19_, exp12_, exp_OCCUSCAN_, etc.
    s = re.sub(r'^exp\d+_', '', s, flags=re.IGNORECASE)
    s = re.sub(r'^exp_[^_]+_', '', s, flags=re.IGNORECASE)
    s = re.sub(r'^exp_', '', s, flags=re.IGNORECASE)
    return s


def load_elisa_binary_csv(path):
    # Read only the first two columns, ignore any extra empty columns, and ignore the header line.
    elisa = pd.read_csv(
        path,
        header=None,
        skiprows=1,
        usecols=[0, 1],
        names=["sample", "elisa_binary"],
        dtype={0: str},
        engine="python",
    )
    elisa["sample"] = elisa["sample"].astype(str).str.strip()
    elisa["elisa_binary"] = pd.to_numeric(elisa["elisa_binary"], errors="coerce")

    # Drop blank / malformed rows
    elisa = elisa.dropna(subset=["sample", "elisa_binary"]).copy()
    elisa = elisa[elisa["sample"].astype(str).str.strip() != ""].copy()

    elisa["elisa_binary"] = elisa["elisa_binary"].astype(int)
    elisa = elisa[elisa["elisa_binary"].isin([0, 1])].copy()

    # Add merge key that strips common run prefixes from PhIP sample names
    elisa["sample_merge"] = elisa["sample"].apply(normalise_sample_for_merge)

    # Sort and keep last if duplicates exist
    elisa = (
        elisa.sort_values(["sample_merge", "sample"])
             .drop_duplicates(subset=["sample_merge"], keep="last")
             .reset_index(drop=True)
    )
    return elisa



def get_samples_to_exclude_from_final_stats(selected_pos_cols, selected_neg_cols):
    # Exclude user-selected comparison / negative controls from final ELISA stats.
    # Keep selected positives in the output by default.
    return set(str(s) for s in selected_neg_cols)

def roc_curve_from_scores(y_true, scores):
    y_true = np.asarray(y_true, dtype=int)
    scores = np.asarray(scores, dtype=float)
    if y_true.size == 0 or len(np.unique(y_true)) < 2:
        return None

    thresholds = np.r_[np.inf, np.unique(scores)[::-1]]
    rows = []
    pos = float((y_true == 1).sum())
    neg = float((y_true == 0).sum())
    for thr in thresholds:
        pred = scores >= thr
        tp = float(((pred == 1) & (y_true == 1)).sum())
        fp = float(((pred == 1) & (y_true == 0)).sum())
        tpr = tp / pos if pos > 0 else 0.0
        fpr = fp / neg if neg > 0 else 0.0
        rows.append((thr, fpr, tpr))

    roc_df = pd.DataFrame(rows, columns=["threshold", "fpr", "tpr"])
    roc_df = roc_df.sort_values(["fpr", "tpr", "threshold"], ascending=[True, True, False]).reset_index(drop=True)
    auc_val = float(np.trapezoid(roc_df["tpr"].values, roc_df["fpr"].values))
    return roc_df, auc_val


def write_roc_outputs(merged_df, prefix, exclude_samples=None):
    usable = merged_df.dropna(subset=["elisa_binary"]).copy()
    if exclude_samples is None:
        exclude_samples = set()
    else:
        exclude_samples = set(str(x) for x in exclude_samples)
    if "sample" in usable.columns and exclude_samples:
        usable = usable[~usable["sample"].astype(str).isin(exclude_samples)].copy()

    metric_map = [
        ("transform1_raw_log10_cpm_sum", "log10(CPM+1) panel sum"),
        ("transform1_raw_log10_cpm_sum_min5", "log10(CPM+1) panel sum, breadth-gated"),
        ("transform1_raw_log10_cpm_sum_x_log_hits", "log10(CPM+1) panel sum × log2(hit count+1)"),
        ("transform1_raw_log10_cpm_sum_x_hits", "log10(CPM+1) panel sum × hit count"),
        ("transform3_minus_mock_mean_floor_log10_sum", "mock mean-corrected log10(CPM+1) panel sum"),
        ("transform4_minus_positive_mock_median_floor_log10_sum", "positive-mock-median-corrected log10(CPM+1) panel sum"),
        ("transform9_minus_mock_median_log10_sum", "mock median-corrected log10(CPM+1) panel sum"),
        ("transform5_minus_mock_q75_floor_log10_sum", "mock Q75-corrected log10(CPM+1) panel sum"),
        ("transform6_minus_mock_q100_floor_log10_sum", "mock Q100-corrected log10(CPM+1) panel sum"),
        ("transform7_minus_q100_plus_median_floor_log10_sum", "mock Q100+median-corrected log10(CPM+1) panel sum"),
        ("transform8_strict_mean3_5sd_only_log10_sum", "mean_3_5sd thresholded log10(CPM+1) panel sum"),
        ("transform_selected_method_continuous_sum", "selected-method continuous panel sum"),
        ("sum_log10_cpm_panel_hits_only_mean_3_5sd", "mean_3_5sd hits-only log10(CPM+1) sum"),
        ("sum_cpm_panel_hits_only_mean_3_5sd", "mean_3_5sd hits-only CPM sum"),
        ("sum_log10_cpm_panel_hits_only_q100_plus_median", "q100_plus_median hits-only log10(CPM+1) sum"),
        ("sum_cpm_panel_hits_only_q100_plus_median", "q100_plus_median hits-only CPM sum"),
        ("transform2_mock_zscore_sum", "mock z-score panel sum"),
    ]

    auc_rows = []
    curve_frames = []
    plotted_any = False
    plt.figure(figsize=(8, 8))
    for metric, label in metric_map:
        if metric not in usable.columns:
            continue
        sub = usable.dropna(subset=[metric, "elisa_binary"]).copy()
        if sub.empty or sub["elisa_binary"].nunique() < 2:
            continue
        roc_res = roc_curve_from_scores(sub["elisa_binary"].astype(int).values, sub[metric].astype(float).values)
        if roc_res is None:
            continue
        roc_df, auc_val = roc_res
        roc_df.insert(0, "metric", metric)
        roc_df.insert(1, "label", label)
        curve_frames.append(roc_df)
        auc_rows.append({"metric": metric, "label": label, "n": len(sub), "auc": auc_val})
        plt.plot(roc_df["fpr"].values, roc_df["tpr"].values, label=f"{label} (AUC={auc_val:.3f})")
        plotted_any = True

    plt.plot([0, 1], [0, 1], linestyle="--", linewidth=1)
    plt.xlabel("False positive rate")
    plt.ylabel("True positive rate")
    plt.title(f"{prefix}: ROC curves vs binary ELISA")
    if plotted_any:
        plt.legend(fontsize=8, loc="lower right")
    plt.tight_layout()
    roc_plot_out = f"{prefix}_roc_curves.png"
    plt.savefig(roc_plot_out, dpi=300)
    plt.close()

    roc_auc_out = f"{prefix}_roc_auc.csv"
    pd.DataFrame(auc_rows, columns=["metric", "label", "n", "auc"]).sort_values("auc", ascending=False).to_csv(roc_auc_out, index=False)

    roc_points_out = f"{prefix}_roc_points.csv"
    if curve_frames:
        pd.concat(curve_frames, ignore_index=True).to_csv(roc_points_out, index=False)
    else:
        pd.DataFrame(columns=["metric", "label", "threshold", "fpr", "tpr"]).to_csv(roc_points_out, index=False)

    return roc_plot_out, roc_auc_out, roc_points_out


def write_orf_profile_plot(matrix_df, sample_cols, prefix, merged_df):
    if "elisa_binary" not in merged_df.columns:
        return None

    merge_key = merged_df[["sample", "elisa_binary"]].dropna().copy()
    if merge_key.empty:
        return None

    ordered_samples = [s for s in merge_key["sample"].astype(str) if s in sample_cols]
    if not ordered_samples:
        return None

    work = matrix_df.copy()
    work["peptide_location"] = work["peptide"].astype(str).map(extract_peptide_location)
    work["location_start"] = work["peptide_location"].map(lambda x: location_sort_key(x)[0])
    work["location_end"] = work["peptide_location"].map(lambda x: location_sort_key(x)[1])

    # Drop rows with invalid or sentinel coordinates before plotting so we do not
    # attempt to allocate a coordinate matrix out to 10**12 when a location field
    # cannot be parsed.
    loc_start_num = pd.to_numeric(work["location_start"], errors="coerce")
    loc_end_num = pd.to_numeric(work["location_end"], errors="coerce")
    valid_coord_mask = (
        loc_start_num.notna() &
        loc_end_num.notna() &
        (loc_start_num >= 0) &
        (loc_end_num >= 0) &
        (loc_start_num < 10**11) &
        (loc_end_num < 10**11)
    )
    invalid_rows = work.loc[~valid_coord_mask, ["peptide", "peptide_location", "location_start", "location_end"]].copy()
    if not invalid_rows.empty:
        invalid_out = f"{prefix}_orf_plot_invalid_coordinates.csv"
        invalid_rows.to_csv(invalid_out, index=False)
        print(f"WARNING: Excluding {len(invalid_rows)} peptide(s) with invalid coordinates from ORF plots: {invalid_out}")
    work = work.loc[valid_coord_mask].copy()
    if work.empty:
        print("WARNING: No valid coordinate-bearing peptides available for ORF profile plot.")
        return None

    work["location_start"] = pd.to_numeric(work["location_start"], errors="coerce").astype(int)
    work["location_end"] = pd.to_numeric(work["location_end"], errors="coerce").astype(int)
    work = work.sort_values(["location_start", "location_end", "peptide_location"]).reset_index(drop=True)

    starts = work["location_start"].astype(int).values
    ends = work["location_end"].astype(int).values
    # Use start + 55 aa where needed, as requested, while respecting larger explicit ends.
    ends = np.maximum(ends, starts + 55)
    max_end = int(np.nanmax(ends)) if len(ends) else 0
    if max_end <= 0:
        print("WARNING: No positive coordinate range available for ORF profile plot.")
        return None

    # Build a continuous coordinate-by-sample hit matrix using log10(CPM+1) values.
    values = np.log10(work[ordered_samples].apply(pd.to_numeric, errors="coerce").fillna(0.0) + 1.0).values
    coord_matrix = np.zeros((len(ordered_samples), max_end + 1), dtype=float)
    hit_count_matrix = np.zeros((len(ordered_samples), max_end + 1), dtype=int)

    for row_i in range(len(work)):
        s = max(0, int(starts[row_i]))
        e = min(max_end, int(ends[row_i]))
        if e < s:
            continue
        segment = slice(s, e + 1)
        row_vals = values[row_i, :]
        for sample_i, v in enumerate(row_vals):
            if v > 0:
                coord_matrix[sample_i, segment] += float(v)
                hit_count_matrix[sample_i, segment] += 1

    # Convert from summed intensity to number of hits per location per sample, as requested.
    hit_plot = hit_count_matrix.astype(float)

    fig, ax = plt.subplots(figsize=(max(16, max_end * 0.03), max(8, len(ordered_samples) * 0.22)))
    im = ax.imshow(hit_plot, aspect='auto', interpolation='nearest', origin='upper')
    ax.set_xlim(0, max_end)
    ax.set_xlabel('Coordinate')
    ax.set_ylabel('Samples')
    ax.set_yticks(range(len(ordered_samples)))
    ax.set_yticklabels(ordered_samples, fontsize=6)
    ax.set_title(f'{prefix}: hits per location per sample (start to start+55 aa / explicit end)')
    cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
    cbar.set_label('Number of peptide hits covering coordinate')
    plt.tight_layout()
    out = f"{prefix}_orf_profile_hits_per_sample.png"
    plt.savefig(out, dpi=300)
    plt.close(fig)

    # Also write the requested group summary over the same continuous x-axis.
    pos_samples = [s for s in merge_key.loc[merge_key["elisa_binary"] == 1, "sample"].astype(str) if s in ordered_samples]
    neg_samples = [s for s in merge_key.loc[merge_key["elisa_binary"] == 0, "sample"].astype(str) if s in ordered_samples]
    pos_idx = [ordered_samples.index(s) for s in pos_samples]
    neg_idx = [ordered_samples.index(s) for s in neg_samples]
    if pos_idx and neg_idx:
        pos_med = np.median(coord_matrix[pos_idx, :], axis=0)
        neg_med = np.median(coord_matrix[neg_idx, :], axis=0)
        x = np.arange(max_end + 1)
        fig2, ax2 = plt.subplots(figsize=(max(16, max_end * 0.03), 5))
        ax2.plot(x, pos_med, label='ELISA positive median log10(CPM+1)')
        ax2.plot(x, neg_med, label='ELISA negative median log10(CPM+1)')
        ax2.set_xlim(0, max_end)
        ax2.set_xlabel('Coordinate')
        ax2.set_ylabel('Median log10(CPM + 1)')
        ax2.set_title(f'{prefix}: continuous coordinate median log10(CPM+1)')
        ax2.legend()
        plt.tight_layout()
        out2 = f"{prefix}_orf_profile_coordinate_median_logcpm.png"
        plt.savefig(out2, dpi=300)
        plt.close(fig2)
    else:
        out2 = None

    profile_csv = f"{prefix}_orf_profile_hits_per_sample.csv"
    pd.DataFrame(hit_plot, index=ordered_samples, columns=list(range(max_end + 1))).to_csv(profile_csv)
    return out, profile_csv


def summarise_elisa_relationships(merged_df, out_csv, exclude_samples=None):
    rows = []
    usable = merged_df.dropna(subset=["elisa_binary"]).copy()
    if exclude_samples is None:
        exclude_samples = set()
    else:
        exclude_samples = set(str(x) for x in exclude_samples)
    if "sample" in usable.columns and exclude_samples:
        usable = usable[~usable["sample"].astype(str).isin(exclude_samples)].copy()
    empty_cols = ["metric", "comparison", "n", "statistic", "p_value", "notes"]
    if usable.empty:
        out_df = pd.DataFrame(columns=empty_cols)
        out_df.to_csv(out_csv, index=False)
        return out_df

    def add_spearman(metric_col, label):
        sub = usable.dropna(subset=[metric_col])
        if sub.empty:
            return
        rho, p = spearmanr(sub[metric_col], sub["elisa_binary"])
        rows.append({"metric": metric_col, "comparison": f"Spearman vs binary ELISA ({label})", "n": len(sub), "statistic": rho, "p_value": p, "notes": label})

    def add_mw(metric_col, label):
        sub = usable.dropna(subset=[metric_col])
        g1 = sub[sub["elisa_binary"] == 1][metric_col]
        g0 = sub[sub["elisa_binary"] == 0][metric_col]
        if len(g1) == 0 or len(g0) == 0:
            return
        stat, p = mannwhitneyu(g1, g0, alternative="two-sided")
        rows.append({"metric": metric_col, "comparison": f"Mann-Whitney 1 vs 0 ({label})", "n": len(sub), "statistic": stat, "p_value": p, "notes": label})

    def add_fisher(call_col, label):
        sub = usable.dropna(subset=[call_col])
        if sub.empty:
            return
        a = int(((sub[call_col] == 1) & (sub["elisa_binary"] == 1)).sum())
        b = int(((sub[call_col] == 1) & (sub["elisa_binary"] == 0)).sum())
        c = int(((sub[call_col] == 0) & (sub["elisa_binary"] == 1)).sum())
        d = int(((sub[call_col] == 0) & (sub["elisa_binary"] == 0)).sum())
        if (a + b + c + d) == 0:
            return
        odds, p = fisher_exact([[a, b], [c, d]])
        rows.append({"metric": call_col, "comparison": f"Fisher exact vs binary ELISA ({label})", "n": len(sub), "statistic": odds, "p_value": p, "notes": f"TP={a}; FP={b}; FN={c}; TN={d}"})

    continuous_metrics = [
        ("transform1_raw_log10_cpm_sum", "transform 1: log10(CPM+1) panel sum"),
        ("transform2_mock_zscore_sum", "transform 2: z scores versus mocks panel sum"),
        ("transform1_raw_log10_cpm_sum_min5", "transform 1: log10(CPM+1) panel sum, zeroed unless >= selected minimum peptides above CPM threshold"),
        ("transform1_raw_log10_cpm_sum_x_log_hits", "transform 1: log10(CPM+1) panel sum multiplied by log2(hit_count+1), with hit_count zeroed unless >= selected minimum peptides above CPM threshold"),
        ("transform1_raw_log10_cpm_sum_x_hits", "transform 1: log10(CPM+1) panel sum multiplied by hit_count, with hit_count zeroed unless >= selected minimum peptides above CPM threshold"),
        ("transform3_minus_mock_mean_floor_log10_sum", "mock mean-corrected summed log10(CPM+1)"),
        ("transform4_minus_positive_mock_median_floor_log10_sum", "positive-mock-median-corrected summed log10(CPM+1)"),
        ("transform9_minus_mock_median_log10_sum", "mock median-corrected summed log10(CPM+1)"),
        ("transform5_minus_mock_q75_floor_log10_sum", "mock Q75-corrected summed log10(CPM+1)"),
        ("transform6_minus_mock_q100_floor_log10_sum", "mock Q100-corrected summed log10(CPM+1)"),
        ("transform7_minus_q100_plus_median_floor_log10_sum", "mock Q100+median-corrected summed log10(CPM+1)"),
        ("transform8_strict_mean3_5sd_only_log10_sum", "mean_3_5sd-thresholded summed log10(CPM+1)"),
        ("transform_selected_method_continuous_sum", "selected-method continuous panel sum"),
        ("sum_log10_cpm_panel_hits_only_mean_3_5sd", "mean_3_5sd hits-only summed log10(CPM+1)"),
        ("sum_cpm_panel_hits_only_mean_3_5sd", "mean_3_5sd hits-only summed CPM"),
        ("sum_log10_cpm_panel_hits_only_q100_plus_median", "q100_plus_median hits-only summed log10(CPM+1)"),
        ("sum_cpm_panel_hits_only_q100_plus_median", "q100_plus_median hits-only summed CPM"),
    ]
    for metric, label in continuous_metrics:
        if metric in usable.columns:
            add_spearman(metric, label)
            add_mw(metric, label)

    binary_metrics = [
        ("call_mean_3_5sd_binary_relaxed", "mean_3_5sd relaxed binary call"),
        ("call_mean_3_5sd_binary_strict", "mean_3_5sd strict binary call"),
        ("call_q100_plus_median_binary_relaxed", "q100_plus_median relaxed binary call"),
        ("call_q100_plus_median_binary_strict", "q100_plus_median strict binary call"),
        ("call_hybrid_binary_relaxed", "hybrid relaxed binary call"),
        ("call_hybrid_binary_strict", "hybrid strict binary call"),
        ("z_mock_binary_relaxed", "z-score relaxed binary call"),
        ("z_mock_binary_strict", "z-score strict binary call"),
    ]
    for metric, label in binary_metrics:
        if metric in usable.columns:
            add_fisher(metric, label)

    out_df = pd.DataFrame(rows, columns=empty_cols)
    out_df.to_csv(out_csv, index=False)
    return out_df


parser = argparse.ArgumentParser(description="PhIP-seq merged analysis and heatmap generation")
parser.add_argument("matrix_file", help="CPM matrix CSV. This can be a full matrix; use --filter-term to subset peptides.")
parser.add_argument("virus_prefix", nargs="?", help="Backward-compatible filter/output term, e.g. 'Africa2 AND |s|'. If --filter-term is not supplied, this is used as the peptide filter.")
parser.add_argument("--filter-term", help="Peptide filter term. Supports AND plus simple OR groups, e.g. 'cchfv AND (|s| OR |m|) AND Africa3'.")
parser.add_argument("--output-prefix", help="Clean output prefix. If omitted, derived from virus_prefix/filter-term or the input filename.")
parser.add_argument("--default", action="store_true", help="Run non-interactively using the standard default preset")
parser.add_argument("--strict", action="store_true", help="Run non-interactively using the strict preset")
parser.add_argument("--common-virus", action="store_true", help="Treat as a common/high-prevalence virus and derive thresholds from mocks")
parser.add_argument("--supervised", action="store_true", help="Run in supervised mode using known positive and negative samples")
parser.add_argument("--unsupervised", action="store_true", help="Run in unsupervised mode without known positive samples")
parser.add_argument("--mock-regex", help="Regex for mock samples")
parser.add_argument("--pos-regex", help="Regex for positive samples; quote alternation patterns such as 'epstein|EBV'")
parser.add_argument("--neg-regex", help="Regex for negative comparison samples; quote alternation patterns such as 'HC|pool'")
parser.add_argument("--pos-select", help="Comma-separated 1-based sample indices (or ranges like 1-3,8) to mark as known positives")
parser.add_argument("--neg-select", help="Comma-separated 1-based sample indices (or ranges like 4-6,9) to mark as known negatives")
parser.add_argument("--heatmap-mode", help="Heatmap mode name or number")
parser.add_argument("--call-mode", help="Primary call to surface in the generic call columns: mean_3_5sd, q100_plus_median, hybrid, both, or all")
parser.add_argument("--include-mocks", action="store_true", help="Force inclusion of mock columns in heatmap display")
parser.add_argument("--mocks", action="store_true", help="Alias for --include-mocks")
parser.add_argument("--exclude-mocks", action="store_true", help="Force exclusion of mock columns in heatmap display")
parser.add_argument("--min-positive-peptides", type=int, help="Minimum positive peptides required to display a sample column")
parser.add_argument("--panel-sizes", help="Comma-separated panel sizes to test")
parser.add_argument("--cpm-thresholds", help="Comma-separated peptide CPM thresholds to test")
parser.add_argument("--panel-choice", type=int, help="Final panel size to use")
parser.add_argument("--cpm-choice", type=float, help="Final peptide CPM threshold to use")
parser.add_argument("--strict-hit-thr", type=int, help="Manual mean_3_5sd hit threshold")
parser.add_argument("--strict-sum-thr", type=float, help="Manual mean_3_5sd hits-only sumCPM threshold")
parser.add_argument("--q100med-hit-thr", type=int, help="Manual Q100+median hit threshold")
parser.add_argument("--q100med-sum-thr", type=float, help="Manual Q100+median sumCPM threshold")
parser.add_argument("--positive-peptides", type=int, help="Use one shared positive-peptide threshold for both calling hit thresholds and heatmap sample suppression")
parser.add_argument("--call-score", choices=["raw_cpm", "mock_subtracted"], default="raw_cpm",
                    help="Score the hybrid call on the raw log10-CPM sum over the panel (raw_cpm, the "
                         "historical behaviour) or on the mock-subtracted transform chosen by "
                         "--heatmap-mode, with breadth counted above the per-peptide mock 3.5SD "
                         "threshold (mock_subtracted). raw_cpm sums background shared by every well, "
                         "so on a large panel the score scale is set by panel size rather than by "
                         "signal and the negative-derived thresholds can land above the highest "
                         "score any sample can reach.")
parser.add_argument("--hybrid-min-negative-hits", type=int, help="Minimum hit count below which the hybrid call is forced negative/non-reactive")
parser.add_argument("--hybrid-min-positive-hits", type=int, help="Minimum hit count required for a hybrid positive/reactive call")
parser.add_argument("--z-clip", type=float, help="Z-score display cap")
parser.add_argument("--output-set", choices=["minimal", "full"], default="minimal", help="minimal writes the core outputs only; full writes all ancillary outputs")
parser.add_argument("--elisa-file", help="Optional CSV file: sample names in column 1 and binary ELISA result (1/0) in column 2; header line is ignored")
parser.add_argument("--collapse-by-location-max", action="store_true", help="Within each peptide location, keep only the maximum value per sample and ignore the other peptides at that location")
parser.add_argument("--no-collapse-by-location-max", action="store_true", help="Do not collapse peptides that share the same location")
parser.add_argument("--collapse-selection-mode", choices=["raw_cpm", "minus_mean_floor", "minus_q75_floor", "minus_q100_floor", "minus_q100_plus_median_floor", "zscore"], help="How to choose the winning peptide within a shared location when collapse is enabled")
parser.add_argument("--collapse-selection-floor", type=float, help="Residual floor to subtract in mock-aware collapse modes")
args = parser.parse_args()

shared_positive_peptides = args.positive_peptides

if args.default and args.strict:
    parser.error("Use either --default or --strict, not both.")
if args.supervised and args.unsupervised:
    parser.error("Use either --supervised or --unsupervised, not both.")
if args.collapse_by_location_max and args.no_collapse_by_location_max:
    parser.error("Use either --collapse-by-location-max or --no-collapse-by-location-max, not both.")

preset = get_preset(use_default=args.default, use_strict=args.strict)

matrix_file = args.matrix_file
# Separate the biological peptide filter from the output prefix.
# Backward compatibility: if --filter-term is omitted, the positional virus_prefix
# still acts as both filter term and output prefix, as in earlier versions.
filter_term = args.filter_term if args.filter_term is not None else args.virus_prefix
prefix_source = args.output_prefix if args.output_prefix is not None else args.virus_prefix
virus_prefix = derive_prefix(matrix_file, prefix_source)

print("\nPhIP-seq merged analysis v1.7.0")
print("This script will:")
print("  1) rank discriminating peptides")
print("  2) tune panel size and peptide CPM threshold")
print("  3) call seropositivity / reactivity")
print("  4) plot a heatmap using the selected panel")
print("")
print(f"Input matrix: {matrix_file}")
print(f"Output prefix: {virus_prefix}")
print(f"Strict mock-based thresholding uses mean + {STRICT_SD_MULTIPLIER} SD")
if preset is not None:
    preset_name = "strict" if args.strict else "default"
    print(f"Preset mode: {preset_name}")
print("")

df = pd.read_csv(matrix_file)

if "peptide" not in df.columns:
    print("ERROR: The input matrix must contain a 'peptide' column.")
    sys.exit(1)

# ------------------------------------------------------------
# Restrict the matrix to peptides matching all requested search terms,
# then collapse peptides sharing the same trailing coordinates after '|'.
# The per-location matrix stores the maximum signal observed for any
# peptide at that location across all sample columns.
# ------------------------------------------------------------
search_term = filter_term if filter_term else ""
search_term = str(search_term).strip()
search_terms = parse_search_terms(search_term)
filtered_matrix_out = None
collapsed_matrix_out = None

if search_term:
    peptide_series = df["peptide"].astype(str)
    mask, parsed_filter_terms = build_filter_mask(peptide_series, search_term)
    matched_n = int(mask.sum())

    if matched_n == 0:
        print(f"ERROR: No peptides matched peptide filter expression: {search_term}")
        print("Parsed filter terms:", parsed_filter_terms)
        sys.exit(1)

    df = df.loc[mask].copy()
    filtered_matrix_out = f"{virus_prefix}_filtered_input_matrix.csv"
    df.to_csv(filtered_matrix_out, index=False)
    print(f"Filtered input matrix using peptide filter expression: {search_term}")
    print(f"Parsed filter terms: {parsed_filter_terms}")
    print(f"Matched peptides retained before location collapse: {matched_n}")
    print(f"Wrote filtered matrix: {filtered_matrix_out}")

sample_cols = [c for c in df.columns if is_true_sample_column(c)]
# Provisional mock detection is needed before optional location collapse because
# mock-aware collapse selection uses mock columns to score candidate peptides.
provisional_mock_regex = args.mock_regex or (preset["mock_regex"] if preset is not None else "mock")
provisional_mock_pat = re.compile(provisional_mock_regex)
mock_cols = [c for c in sample_cols if provisional_mock_pat.search(c)]
collapsed_matrix_out = None
if sample_cols:
    if args.collapse_by_location_max:
        collapse_by_location_max = True
    elif args.no_collapse_by_location_max:
        collapse_by_location_max = False
    else:
        # Default OFF. Collapsing replaces peptide names with bare position
        # labels, which carry no ORF or coordinates, so the selected panel can no
        # longer be used for genome plots. It also discards real peptides where a
        # library tiles one region from many source entries.
        collapse_by_location_max = yes_no_prompt(
            "If multiple peptides share the same location, keep only the maximum value per sample at that location and ignore the others?",
            default="n"
        )

    if collapse_by_location_max:
        if args.collapse_selection_mode is not None:
            collapse_selection_mode = str(args.collapse_selection_mode).strip().lower()
        else:
            collapse_selection_mode = prompt_with_default(
                "Collapse selection mode (raw_cpm, minus_mean_floor, minus_q75_floor, minus_q100_floor, minus_q100_plus_median_floor, zscore)",
                "raw_cpm"
            ).strip().lower()
        if collapse_selection_mode not in {"raw_cpm", "minus_mean_floor", "minus_q75_floor", "minus_q100_floor", "minus_q100_plus_median_floor", "zscore"}:
            print(f"WARNING: unrecognised collapse selection mode '{collapse_selection_mode}', defaulting to raw_cpm.")
            collapse_selection_mode = "raw_cpm"

        if args.collapse_selection_floor is not None:
            collapse_selection_floor = float(args.collapse_selection_floor)
        elif collapse_selection_mode == "raw_cpm":
            collapse_selection_floor = 0.0
        else:
            collapse_selection_floor = float(prompt_with_default(
                "Residual floor to subtract in mock-aware collapse selection",
                "5"
            ))

        collapsed_df = collapse_matrix_to_locations(
            df,
            sample_cols,
            mock_cols=mock_cols,
            selection_mode=collapse_selection_mode,
            selection_floor=collapse_selection_floor,
        )
        collapsed_matrix_out = f"{virus_prefix}_collapsed_by_location_matrix.csv"
        collapsed_df.to_csv(collapsed_matrix_out, index=False)
        print(f"Collapsed to unique peptide locations using selection mode '{collapse_selection_mode}': {len(collapsed_df)}")
        if collapse_selection_mode != "raw_cpm":
            print(f"Mock-aware collapse residual floor: {collapse_selection_floor}")
        print(f"Wrote collapsed location matrix: {collapsed_matrix_out}")
        df = collapsed_df.copy()
    else:
        print("Keeping all matched peptides, including peptides that share the same location.")

sample_cols = [c for c in df.columns if is_true_sample_column(c)]
if not sample_cols:
    print("ERROR: No sample columns found in the matrix.")
    sys.exit(1)

has_true_negative_controls = yes_no_or_preset(
    "Do you have true negative controls on your run?",
    default="y",
    preset_value=(None if preset is None else (not preset["common_virus_mode"]))
)
common_virus_mode = (not has_true_negative_controls)

mock_regex = args.mock_regex or prompt_or_preset(
    "Enter a regex to identify mock samples (for example 'mock')",
    "mock",
    preset_value=(preset["mock_regex"] if preset else None)
)

mock_pat = re.compile(mock_regex)
mock_cols = [c for c in sample_cols if mock_pat.search(c)]
nonmock_cols = [c for c in sample_cols if c not in mock_cols]

if not mock_cols:
    print("ERROR: No mock samples matched.")
    print("Available sample columns:")
    for s in sample_cols:
        print(" ", s)
    sys.exit(1)

print_indexed_samples(nonmock_cols, "Available NON-MOCK samples (for supervised selection)")

supervised_mode = bool(args.supervised)
if args.unsupervised:
    supervised_mode = False
elif (not args.supervised) and (args.pos_regex or args.neg_regex or args.pos_select or args.neg_select):
    supervised_mode = True
elif not args.supervised:
    supervised_mode = yes_no_prompt(
        "Run in supervised mode using known positive and negative samples?",
        default="n"
    )

pos_regex = args.pos_regex
neg_regex = args.neg_regex
selected_pos_cols = []
selected_neg_cols = []

# Treat blank preset / CLI values as missing so supervised mode can still prompt.
if pos_regex is not None and str(pos_regex).strip() == "":
    pos_regex = None
if neg_regex is not None and str(neg_regex).strip() == "":
    neg_regex = None

if supervised_mode and (args.pos_select or args.neg_select or (preset is None and pos_regex is None and neg_regex is None)):
    if args.neg_select:
        neg_select_text = args.neg_select
    elif preset is not None:
        neg_select_text = ""
    else:
        neg_select_text = input(
            "Enter known NEGATIVE sample indices (comma-separated; ranges allowed, e.g. 1-3,8).\n"
        ).strip()
    if args.pos_select:
        pos_select_text = args.pos_select
    elif preset is not None:
        pos_select_text = ""
    else:
        pos_select_text = input(
            "Enter known POSITIVE sample indices (comma-separated; ranges allowed, e.g. 4-6,9).\n"
        ).strip()
    try:
        selected_neg_cols = parse_index_selection(neg_select_text, nonmock_cols)
        selected_pos_cols = parse_index_selection(pos_select_text, nonmock_cols)
    except ValueError as e:
        print(f"ERROR: {e}")
        sys.exit(1)
    supervised_mode = True
elif not supervised_mode:
    pos_regex = ""

if supervised_mode:
    if pos_regex is None and not selected_pos_cols:
        pos_regex = input(
            "Enter a regex to identify POSITIVE samples, or press Enter to use the index-based selections above.\n"
            "Example: 'MHRA' or 'AV|MHRA'\n"
        ).strip()
    if neg_regex is None and not selected_neg_cols:
        neg_regex = input(
            "Enter a regex to identify NEGATIVE comparison samples, or press Enter to use the index-based selections above.\n"
            "Example: 'HCV_pool|HC17|KOVAG|PRIVIGEN|LDV_IgG|Ig_depleted_serum'\n"
        ).strip()
else:
    pos_regex = ""
    if args.neg_regex is not None:
        neg_regex = args.neg_regex
    elif preset is not None:
        neg_regex = preset["neg_regex"]
    else:
        neg_regex = input(
            "Optional: enter a regex to identify comparison / likely-negative samples for labelling only.\n"
            "Press Enter to skip this.\n"
        ).strip()

unsupervised_mode = not supervised_mode
if unsupervised_mode:
    print("\nRunning in UNSUPERVISED mock-based mode.")
    print("Peptides will be ranked by non-mock enrichment above mocks.")
else:
    print("\nRunning in SUPERVISED mode.")
    if common_virus_mode:
        print("Common-virus mode selected.")
        print("Final thresholds will still be derived from MOCKS, not from named biological negatives.")

pos_pat = re.compile(pos_regex) if pos_regex else None
neg_pat = re.compile(neg_regex) if neg_regex else None

regex_pos_cols = [c for c in sample_cols if pos_pat and pos_pat.search(c)]
regex_neg_cols = [c for c in sample_cols if neg_pat and neg_pat.search(c)]
pos_cols = list(dict.fromkeys(selected_pos_cols + regex_pos_cols))
neg_cols = list(dict.fromkeys(selected_neg_cols + regex_neg_cols))

# keep groups disjoint and exclude mocks
pos_cols = [c for c in pos_cols if c in nonmock_cols]
neg_cols = [c for c in neg_cols if c in nonmock_cols and c not in pos_cols]

if supervised_mode and not pos_cols:
    print("ERROR: Supervised mode requires at least one known positive sample.")
    print("Available non-mock sample columns:")
    for i, s in enumerate(nonmock_cols, start=1):
        print(f"  {i:>3}. {s}")
    sys.exit(1)
if supervised_mode and (not common_virus_mode) and (not neg_cols):
    print("ERROR: Supervised mode requires at least one known negative sample unless common-virus mode is being used.")
    print("Available non-mock sample columns:")
    for i, s in enumerate(nonmock_cols, start=1):
        print(f"  {i:>3}. {s}")
    sys.exit(1)

if args.panel_sizes:
    candidate_sizes = [int(x.strip()) for x in args.panel_sizes.split(",") if x.strip()]
elif preset is not None:
    candidate_sizes = [int(x) for x in preset["candidate_sizes"]]
else:
    candidate_sizes_in = prompt_with_default(
        "Enter peptide panel sizes to test, comma-separated",
        "5,10,15,20,25,30,40,50"
    )
    candidate_sizes = [int(x.strip()) for x in candidate_sizes_in.split(",") if x.strip()]

if args.cpm_thresholds:
    candidate_cpms = [float(x.strip()) for x in args.cpm_thresholds.split(",") if x.strip()]
elif preset is not None:
    candidate_cpms = [float(x) for x in preset["candidate_cpms"]]
else:
    candidate_cpms_in = prompt_with_default(
        "Enter peptide CPM cut-offs to test, comma-separated\n"
        "This is the CPM needed for a peptide to count as present in a sample",
        "10,15,20,30"
    )
    candidate_cpms = [float(x.strip()) for x in candidate_cpms_in.split(",") if x.strip()]

# v1.7: resolve the actual CPM threshold before peptide ranking.
# This threshold is used for mock-aware ranking, panel scoring, calls and R plotting.
if args.cpm_choice is not None:
    cpm_choice = float(args.cpm_choice)
else:
    default_cpm_for_ranking = candidate_cpms[0] if candidate_cpms else 10.0
    cpm_choice = float(prompt_with_default(
        "What CPM cut-off should be used for mock-aware ranking and peptide presence?",
        str(int(default_cpm_for_ranking) if float(default_cpm_for_ranking).is_integer() else default_cpm_for_ranking)
    ))

# In v1.7, panel tuning is performed at the actual chosen CPM threshold only.
# This prevents peptide ranking/selection being optimised at a threshold that the user
# subsequently ignores.
candidate_cpms = [cpm_choice]

print("\nMatched mock samples:")
for s in mock_cols:
    print(" ", s)

if pos_cols:
    print("\nMatched positive samples:")
    for s in pos_cols:
        print(" ", s)

if neg_cols:
    print("\nMatched comparison / negative-labelled samples:")
    for s in neg_cols:
        print(" ", s)

other_nonmock_cols = [c for c in nonmock_cols if c not in pos_cols and c not in neg_cols]
if other_nonmock_cols:
    print("\nOther non-mock / unlabelled samples:")
    for s in other_nonmock_cols:
        print(" ", s)

# rank peptides using the actual selected CPM threshold and selected signal definition
# -----------------------------------------------------------------------------
ranking_method = args.heatmap_mode or (preset["heatmap_mode"] if preset else "mock_q100_plus_median")
if str(ranking_method).isdigit():
    ranking_method = {
        "1": "raw_cpm",
        "2": "zscore",
        "3": "minus_mean",
        "4": "positive_mock_median",
        "5": "mock_q75",
        "6": "mock_q100",
        "7": "mock_q100_plus_median",
        "8": "mean_3_5sd",
    }.get(str(ranking_method), "mock_q100_plus_median")
ranking_method = str(ranking_method).strip().lower()
if ranking_method == "all":
    ranking_method = "mock_q100_plus_median"

print(f"\nRanking peptides using method '{ranking_method}' and CPM floor {cpm_choice:g}")


def _sample_signal_for_ranking(row, sample_name, method, cpm_floor):
    v = pd.to_numeric(row.get(sample_name, np.nan), errors="coerce")
    if pd.isna(v):
        return 0.0
    v = float(v)
    stats = get_linked_mock_stats_from_row(pd.Series(row), sample_name, mock_cols)
    cpm_floor = float(cpm_floor)

    # All methods apply the user's chosen CPM floor for ranking so that panel
    # selection and downstream presence calls are aligned.
    if method == "raw_cpm":
        return v if v >= cpm_floor else 0.0

    if method == "zscore":
        if v < cpm_floor:
            return 0.0
        sd = stats["sd"]
        if sd <= 0:
            return 0.0
        return max(0.0, (v - stats["mean"]) / sd)

    if method in {"minus_mean", "mock_mean"}:
        baseline = stats["mean"]
        return max(0.0, v - baseline - cpm_floor)

    if method in {"positive_mock_median", "mock_positive_median"}:
        baseline = stats["positive_median"]
        return max(0.0, v - baseline - cpm_floor)

    if method == "mock_q75":
        baseline = stats["q75"]
        return max(0.0, v - baseline - cpm_floor)

    if method == "mock_q100":
        baseline = stats["q100"]
        return max(0.0, v - baseline - cpm_floor)

    if method == "mean_3_5sd":
        thr = max(cpm_floor, stats["mean"] + STRICT_SD_MULTIPLIER * stats["sd"])
        return v if v > thr else 0.0

    # Default and recommended mock-aware ranking: signal remaining above the
    # maximum mock value, mock median and the selected CPM floor.
    baseline = stats["q100"] + stats["median"]
    return max(0.0, v - baseline - cpm_floor)


rank_signal_cols = []
for sname in sample_cols:
    col = f"__rank_signal__{sname}"
    rank_signal_cols.append(col)
    df[col] = df.apply(lambda row, ss=sname: _sample_signal_for_ranking(row, ss, ranking_method, cpm_choice), axis=1)

nonmock_signal_cols = [f"__rank_signal__{s}" for s in nonmock_cols]
mock_signal_cols = [f"__rank_signal__{s}" for s in mock_cols]
pos_signal_cols = [f"__rank_signal__{s}" for s in pos_cols]
neg_signal_cols = [f"__rank_signal__{s}" for s in neg_cols]

if unsupervised_mode:
    df["nonmock_rank_signal_mean"] = df[nonmock_signal_cols].mean(axis=1)
    df["nonmock_rank_signal_max"] = df[nonmock_signal_cols].max(axis=1)
    df["mock_rank_signal_mean"] = df[mock_signal_cols].mean(axis=1)
    df["mock_rank_signal_max"] = df[mock_signal_cols].max(axis=1)
    above_col = "n_nonmock_positive_rank_signal"
    df[above_col] = df[nonmock_signal_cols].gt(0).sum(axis=1)
    df["nonmock_minus_mock_rank_signal"] = df["nonmock_rank_signal_mean"] - df["mock_rank_signal_mean"]
    df["nonmock_fc_over_mock_rank_signal_plus1"] = (df["nonmock_rank_signal_mean"] + 1.0) / (df["mock_rank_signal_mean"] + 1.0)
    df["score"] = (
        (df[above_col] + 1.0)
        * df["nonmock_fc_over_mock_rank_signal_plus1"]
        * (df["nonmock_minus_mock_rank_signal"] + 1.0).clip(lower=0.0)
    )
    rank_cols = [
        "peptide",
        "nonmock_rank_signal_mean",
        "nonmock_rank_signal_max",
        "mock_rank_signal_mean",
        "mock_rank_signal_max",
        above_col,
        "nonmock_minus_mock_rank_signal",
        "nonmock_fc_over_mock_rank_signal_plus1",
        "score",
    ]
else:
    df["pos_rank_signal_mean"] = df[pos_signal_cols].mean(axis=1)
    df["mock_rank_signal_mean"] = df[mock_signal_cols].mean(axis=1)
    df["mock_rank_signal_max"] = df[mock_signal_cols].max(axis=1)

    if neg_cols:
        df["neg_rank_signal_mean"] = df[neg_signal_cols].mean(axis=1)
        df["positive_minus_negative_rank_signal"] = df["pos_rank_signal_mean"] - df["neg_rank_signal_mean"]
        df["positive_minus_mock_rank_signal"] = df["pos_rank_signal_mean"] - df["mock_rank_signal_mean"]
        df["positive_fc_over_negative_rank_signal_plus1"] = (df["pos_rank_signal_mean"] + 1.0) / (df["neg_rank_signal_mean"] + 1.0)
        df["positive_fc_over_mock_rank_signal_plus1"] = (df["pos_rank_signal_mean"] + 1.0) / (df["mock_rank_signal_mean"] + 1.0)
        df["score"] = (
            df["positive_fc_over_negative_rank_signal_plus1"]
            * df["positive_fc_over_mock_rank_signal_plus1"]
            * (df["positive_minus_mock_rank_signal"] + 1.0).clip(lower=0.0)
        )
        rank_cols = [
            "peptide",
            "pos_rank_signal_mean",
            "neg_rank_signal_mean",
            "mock_rank_signal_mean",
            "mock_rank_signal_max",
            "positive_minus_negative_rank_signal",
            "positive_minus_mock_rank_signal",
            "positive_fc_over_negative_rank_signal_plus1",
            "positive_fc_over_mock_rank_signal_plus1",
            "score",
        ]
    else:
        df["positive_minus_mock_rank_signal"] = df["pos_rank_signal_mean"] - df["mock_rank_signal_mean"]
        df["positive_fc_over_mock_rank_signal_plus1"] = (df["pos_rank_signal_mean"] + 1.0) / (df["mock_rank_signal_mean"] + 1.0)
        df["score"] = (
            df["positive_fc_over_mock_rank_signal_plus1"]
            * (df["positive_minus_mock_rank_signal"] + 1.0).clip(lower=0.0)
        )
        rank_cols = [
            "peptide",
            "pos_rank_signal_mean",
            "mock_rank_signal_mean",
            "mock_rank_signal_max",
            "positive_minus_mock_rank_signal",
            "positive_fc_over_mock_rank_signal_plus1",
            "score",
        ]

# Drop rows with no positive ranking signal across non-mock samples. This is the
# explicit v1.7 behaviour requested: peptides are ranked after mock/background
# removal and after applying the user's CPM threshold.
pre_rank_n = len(df)
df = df[df[nonmock_signal_cols].gt(0).any(axis=1)].copy()
post_rank_n = len(df)
print(f"Peptides retained for ranking after mock-aware CPM filtering: {post_rank_n}/{pre_rank_n}")
if post_rank_n == 0:
    print("ERROR: No peptides retained after mock-aware CPM filtering. Consider lowering CPM threshold or using raw_cpm ranking.")
    sys.exit(1)

ranked_extra_cols = [c for c in ["ambiguous_score", "source_peptide_count", "orig_peptides"] if c in df.columns and c not in rank_cols]
ranked_df = df[rank_cols + ranked_extra_cols].sort_values(
    by=["score"],
    ascending=[False]
).reset_index(drop=True)

ranked_out = f"{virus_prefix}_ranked_peptides.csv"
ranked_df.to_csv(ranked_out, index=False)
print(f"\nWrote ranked peptide table: {ranked_out}")

lookup_cols = list(dict.fromkeys(sample_cols + [c for c in df.columns if "__mock__" in str(c)]))
matrix_lookup = df.set_index("peptide")[lookup_cols].to_dict(orient="index")
ranked_peptides = ranked_df["peptide"].tolist()


def score_panel(top_n, peptide_presence_cpm):
    panel = ranked_peptides[:top_n]
    if not panel:
        return None

    per_sample = compute_per_sample_metrics(panel, peptide_presence_cpm, "mean_3_5sd")

    mock_hits = [r["n_above_mock_mean_3_5sd_panel"] for r in per_sample if r["group"] == "mock"]
    mock_sum = [r["sum_cpm_panel"] for r in per_sample if r["group"] == "mock"]

    if unsupervised_mode:
        comp_hits = [r["n_above_mock_mean_3_5sd_panel"] for r in per_sample if r["group"] != "mock"]
        comp_sum = [r["sum_cpm_panel"] for r in per_sample if r["group"] != "mock"]
        comp_label = "nonmock"
    else:
        pos_hits = [r["n_above_mock_mean_3_5sd_panel"] for r in per_sample if r["group"] == "positive"]
        pos_sum = [r["sum_cpm_panel"] for r in per_sample if r["group"] == "positive"]

        if common_virus_mode:
            comp_hits = mock_hits
            comp_sum = mock_sum
            comp_label = "mock"
        else:
            comp_hits = [r["n_above_mock_mean_3_5sd_panel"] for r in per_sample if r["group"] == "negative"]
            comp_sum = [r["sum_cpm_panel"] for r in per_sample if r["group"] == "negative"]
            comp_label = "negative"

    if unsupervised_mode:
        if not comp_hits or not comp_sum or not mock_hits or not mock_sum:
            return None
        count_sep = (safe_median(comp_hits) - safe_median(mock_hits)) / (safe_iqr(mock_hits) + 1.0)
        cpm_sep = (safe_median(comp_sum) - safe_median(mock_sum)) / (safe_iqr(mock_sum) + 1.0)
        total_sep = count_sep + cpm_sep
        return {
            "panel_size": len(panel),
            "peptide_presence_cpm": peptide_presence_cpm,
            "panel_peptides": panel,
            "per_sample": per_sample,
            "comparison_label": comp_label,
            "signal_median_hits": safe_median(comp_hits),
            "signal_median_sum_cpm": safe_median(comp_sum),
            "mock_median_hits": safe_median(mock_hits),
            "mock_median_sum_cpm": safe_median(mock_sum),
            "count_separation": count_sep,
            "cpm_separation": cpm_sep,
            "mock_penalty": 0.0,
            "total_separation": total_sep,
        }
    else:
        if not pos_hits or not pos_sum or not comp_hits or not comp_sum:
            return None

        count_sep = (safe_median(pos_hits) - safe_median(comp_hits)) / (safe_iqr(comp_hits) + 1.0)
        cpm_sep = (safe_median(pos_sum) - safe_median(comp_sum)) / (safe_iqr(comp_sum) + 1.0)
        mock_penalty = (
            safe_median(mock_sum) / (safe_median(pos_sum) + 1.0)
            + safe_median(mock_hits) / (safe_median(pos_hits) + 1.0)
        )
        total_sep = count_sep + cpm_sep - mock_penalty

        return {
            "panel_size": len(panel),
            "peptide_presence_cpm": peptide_presence_cpm,
            "panel_peptides": panel,
            "per_sample": per_sample,
            "comparison_label": comp_label,
            "pos_median_hits": safe_median(pos_hits),
            "comp_median_hits": safe_median(comp_hits),
            "mock_median_hits": safe_median(mock_hits),
            "pos_median_sum_cpm": safe_median(pos_sum),
            "comp_median_sum_cpm": safe_median(comp_sum),
            "mock_median_sum_cpm": safe_median(mock_sum),
            "count_separation": count_sep,
            "cpm_separation": cpm_sep,
            "mock_penalty": mock_penalty,
            "total_separation": total_sep,
        }


tuning_results = []
for n in candidate_sizes:
    for cpm in candidate_cpms:
        res = score_panel(n, cpm)
        if res is not None:
            tuning_results.append(res)

if not tuning_results:
    print("ERROR: No valid tuning results could be generated.")
    sys.exit(1)

tuning_results.sort(key=lambda x: x["total_separation"], reverse=True)
best = tuning_results[0]

tuning_out = f"{virus_prefix}_panel_tuning.csv"
with open(tuning_out, "w", newline="") as f:
    w = csv.writer(f)
    if unsupervised_mode:
        w.writerow([
            "panel_size",
            "peptide_presence_cpm",
            "comparison_label",
            "signal_median_hits",
            "mock_median_hits",
            "signal_median_sum_cpm",
            "mock_median_sum_cpm",
            "count_separation",
            "cpm_separation",
            "total_separation",
        ])
        for r in sorted(tuning_results, key=lambda x: (x["panel_size"], x["peptide_presence_cpm"])):
            w.writerow([
                r["panel_size"],
                r["peptide_presence_cpm"],
                r["comparison_label"],
                r["signal_median_hits"],
                r["mock_median_hits"],
                round(r["signal_median_sum_cpm"], 4),
                round(r["mock_median_sum_cpm"], 4),
                round(r["count_separation"], 4),
                round(r["cpm_separation"], 4),
                round(r["total_separation"], 4),
            ])
    else:
        w.writerow([
            "panel_size",
            "peptide_presence_cpm",
            "comparison_label",
            "pos_median_hits",
            "comparison_median_hits",
            "mock_median_hits",
            "pos_median_sum_cpm",
            "comparison_median_sum_cpm",
            "mock_median_sum_cpm",
            "count_separation",
            "cpm_separation",
            "mock_penalty",
            "total_separation",
        ])
        for r in sorted(tuning_results, key=lambda x: (x["panel_size"], x["peptide_presence_cpm"])):
            w.writerow([
                r["panel_size"],
                r["peptide_presence_cpm"],
                r["comparison_label"],
                r["pos_median_hits"],
                r["comp_median_hits"],
                r["mock_median_hits"],
                round(r["pos_median_sum_cpm"], 4),
                round(r["comp_median_sum_cpm"], 4),
                round(r["mock_median_sum_cpm"], 4),
                round(r["count_separation"], 4),
                round(r["cpm_separation"], 4),
                round(r["mock_penalty"], 4),
                round(r["total_separation"], 4),
            ])

print("\nTuning results (higher score is better):")
for r in sorted(tuning_results, key=lambda x: (x["panel_size"], x["peptide_presence_cpm"])):
    if unsupervised_mode:
        print(
            f"  panel={r['panel_size']:>3}  "
            f"cpm={int(r['peptide_presence_cpm']):>2}  "
            f"hits(nonmock/mock)={r['signal_median_hits']}/{r['mock_median_hits']}  "
            f"sumCPM(nonmock/mock)={round(r['signal_median_sum_cpm'],2)}/{round(r['mock_median_sum_cpm'],2)}  "
            f"score={round(r['total_separation'],3)}"
        )
    else:
        print(
            f"  panel={r['panel_size']:>3}  "
            f"cpm={int(r['peptide_presence_cpm']):>2}  "
            f"hits(pos/{r['comparison_label']}/mock)="
            f"{r['pos_median_hits']}/{r['comp_median_hits']}/{r['mock_median_hits']}  "
            f"sumCPM(pos/{r['comparison_label']}/mock)="
            f"{round(r['pos_median_sum_cpm'],2)}/{round(r['comp_median_sum_cpm'],2)}/{round(r['mock_median_sum_cpm'],2)}  "
            f"score={round(r['total_separation'],3)}"
        )

tuned_panel_size = int(best["panel_size"])
recommended_panel_size = len(ranked_peptides)

print(f"\nTuned panel size suggestion: {tuned_panel_size} peptides")
print(f"Defaulting to tuned panel size: {tuned_panel_size} peptides (full matched/ranked panel: {recommended_panel_size})")
tuned_cpm_choice = float(best["peptide_presence_cpm"])
tuned_cpm_print = int(tuned_cpm_choice) if float(tuned_cpm_choice).is_integer() else tuned_cpm_choice
print(f"Panel tuning was performed at selected CPM cut-off: {cpm_choice:g}")

if args.panel_choice is not None:
    panel_choice = int(args.panel_choice)
else:
    panel_choice = int(prompt_with_default(
        "How many top peptides would you like to use?",
        str(tuned_panel_size)
    ))

# Write the actual user-selected settings immediately after the final choices are resolved.
# Downstream plotting should use these values, not the tuning recommendation.
settings_out = f"{virus_prefix}_run_settings.csv"
with open(settings_out, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=[
        "output_prefix",
        "filter_term",
        "panel_choice",
        "cpm_choice",
        "tuned_panel_size",
        "tuned_cpm_choice",
        "heatmap_mode",
        "ranking_method",
        "selected_continuous_metric",
    ])
    writer.writeheader()
    writer.writerow({
        "output_prefix": virus_prefix,
        "filter_term": search_term,
        "panel_choice": panel_choice,
        "cpm_choice": cpm_choice,
        "tuned_panel_size": tuned_panel_size,
        "tuned_cpm_choice": tuned_cpm_choice,
        "heatmap_mode": args.heatmap_mode or (preset["heatmap_mode"] if preset else "mock_q100_plus_median"),
        "ranking_method": ranking_method,
        "selected_continuous_metric": selected_continuous_signal_metric_name(ranking_method),
    })
print(f"Wrote run settings: {settings_out}")

final = score_panel(panel_choice, cpm_choice)
if final is None:
    print("ERROR: Could not compute final panel score.")
    sys.exit(1)

# Resolve the shared peptide-count threshold used for breadth-gated metrics,
# heatmap sample suppression, and the default hybrid positive hit threshold.
if shared_positive_peptides is not None:
    derived_metric_min_hits = int(shared_positive_peptides)
elif args.min_positive_peptides is not None:
    derived_metric_min_hits = int(args.min_positive_peptides)
elif preset is not None:
    derived_metric_min_hits = int(preset["min_positive_peptides_for_display"])
else:
    derived_metric_min_hits = int(prompt_with_default(
        "Minimum number of positive peptides / hits to use for breadth-gated metrics and default hybrid positive calls",
        "5"
    ))

# ------------------------------------------------------------
# build initial report table and compute independent qualitative calls
# for mean_3_5sd, q100_plus_median, and hybrid.
# ------------------------------------------------------------
mean_metrics = compute_per_sample_metrics(final["panel_peptides"], cpm_choice, "mean_3_5sd")
q100med_metrics = compute_per_sample_metrics(final["panel_peptides"], cpm_choice, "q100_plus_median")
mean_hit_col = "n_above_mock_mean_3_5sd_panel"
q100med_hit_col = "n_above_mock_q100_plus_median_panel"
mean_sum_hits_col = "sum_cpm_panel_hits_only_mean_3_5sd"
q100med_sum_hits_col = "sum_cpm_panel_hits_only_q100_plus_median"
mean_by_sample = {r["sample"]: r for r in mean_metrics}
q100med_by_sample = {r["sample"]: r for r in q100med_metrics}

report_rows = []
for s in sample_cols:
    base_group = mean_by_sample[s]["group"] if s in mean_by_sample else q100med_by_sample[s]["group"]
    row = {
        "sample": s,
        "group": base_group,
        "n_panel_peptides": len(final["panel_peptides"]),
        "sum_cpm_panel": mean_by_sample[s]["sum_cpm_panel"],
        "n_above_mock_mean_3_5sd_panel": mean_by_sample[s][mean_hit_col],
        "sum_cpm_panel_hits_only_mean_3_5sd": mean_by_sample[s][mean_sum_hits_col],
        "n_above_mock_q100_plus_median_panel": q100med_by_sample[s][q100med_hit_col],
        "sum_cpm_panel_hits_only_q100_plus_median": q100med_by_sample[s][q100med_sum_hits_col],
    }
    report_rows.append(row)

report_out = f"{virus_prefix}_seropositivity_report.csv"
report_df = pd.DataFrame(report_rows)
report_df = report_df.sort_values(["sample"]).reset_index(drop=True)

# Add thresholded CPM and mock-based z-score metrics using the selected panel
z_metrics = compute_mock_zscore_metrics(df, final["panel_peptides"], sample_cols, mock_cols, cpm_choice)
for colname, mapping in z_metrics.items():
    report_df[colname] = report_df["sample"].map(mapping)

# Add hit-only summed metrics aligned to the hit definitions
hits_only_mean = compute_hits_only_metrics(df, final["panel_peptides"], sample_cols, mock_cols, cpm_choice, "mean_3_5sd")
hits_only_q100 = compute_hits_only_metrics(df, final["panel_peptides"], sample_cols, mock_cols, cpm_choice, "q100_plus_median")
for metric_dict in [hits_only_mean, hits_only_q100]:
    for colname, mapping in metric_dict.items():
        if colname not in report_df.columns:
            report_df[colname] = report_df["sample"].map(mapping)

# Add panel summaries for all heatmap-style transforms (1-8)
transform_metrics = compute_heatmap_transform_panel_metrics(
    df, final["panel_peptides"], sample_cols, mock_cols, cpm_choice, None
)
for colname, mapping in transform_metrics.items():
    report_df[colname] = report_df["sample"].map(mapping)

selected_continuous_metric = selected_continuous_signal_metric_name(ranking_method)
if selected_continuous_metric in report_df.columns:
    report_df["transform_selected_method_continuous_sum"] = report_df[selected_continuous_metric]
else:
    report_df["transform_selected_method_continuous_sum"] = np.nan
report_df["selected_continuous_metric_source"] = selected_continuous_metric

# Derived logCPM metrics based on the number of selected panel peptides above the chosen CPM threshold.
panel_df_for_hits = df[df["peptide"].isin(final["panel_peptides"])].copy()
if not panel_df_for_hits.empty:
    hit_counts = {}
    for sample in sample_cols:
        vals = pd.to_numeric(panel_df_for_hits[sample], errors="coerce").fillna(0)
        hit_counts[sample] = int((vals >= float(cpm_choice)).sum())
    report_df["panel_peptide_hit_count"] = report_df["sample"].map(hit_counts).fillna(0).astype(int)
else:
    report_df["panel_peptide_hit_count"] = 0

effective_hit_counts = np.where(
    report_df["panel_peptide_hit_count"] >= derived_metric_min_hits,
    report_df["panel_peptide_hit_count"].astype(float),
    0.0
)
report_df["transform1_raw_log10_cpm_sum_min5"] = np.where(
    report_df["panel_peptide_hit_count"] >= derived_metric_min_hits,
    report_df["transform1_raw_log10_cpm_sum"],
    0.0
)
report_df["transform1_raw_log10_cpm_sum_x_log_hits"] = (
    report_df["transform1_raw_log10_cpm_sum"] *
    np.log2(effective_hit_counts + 1.0)
)
report_df["transform1_raw_log10_cpm_sum_x_hits"] = (
    report_df["transform1_raw_log10_cpm_sum"] *
    effective_hit_counts
)

# Independent mean_3_5sd and q100_plus_median qualitative calls.
use_mocks_for_thresholds = bool(common_virus_mode or unsupervised_mode)

mean_sum_log_hits_col = "sum_log10_cpm_panel_hits_only_mean_3_5sd"
q100_sum_log_hits_col = "sum_log10_cpm_panel_hits_only_q100_plus_median"

mean_thresholds = derive_call_thresholds(mean_metrics, mean_hit_col, mean_sum_log_hits_col, use_mocks_for_thresholds)
q100_thresholds = derive_call_thresholds(q100med_metrics, q100med_hit_col, q100med_sum_hits_col, use_mocks_for_thresholds)

mean_hit_thr = int(args.strict_hit_thr) if args.strict_hit_thr is not None else int(mean_thresholds["hit_threshold_default"])
mean_sum_thr = float(args.strict_sum_thr) if args.strict_sum_thr is not None else float(mean_thresholds["sum_threshold_default"])
q100_hit_thr = int(args.q100med_hit_thr) if args.q100med_hit_thr is not None else int(q100_thresholds["hit_threshold_default"])
q100_sum_thr = float(args.q100med_sum_thr) if args.q100med_sum_thr is not None else float(q100_thresholds["sum_threshold_default"])

mean_called = apply_calls_to_metrics(mean_metrics, mean_hit_col, mean_sum_log_hits_col, mean_hit_thr, mean_sum_thr, unsupervised_mode)
q100_called = apply_calls_to_metrics(q100med_metrics, q100med_hit_col, q100med_sum_hits_col, q100_hit_thr, q100_sum_thr, unsupervised_mode)

mean_call_by_sample = {r["sample"]: r["call"] for r in mean_called}
q100_call_by_sample = {r["sample"]: r["call"] for r in q100_called}

report_df["call_mean_3_5sd"] = report_df["sample"].map(mean_call_by_sample)
report_df["call_q100_plus_median"] = report_df["sample"].map(q100_call_by_sample)
report_df["mean_3_5sd_hit_threshold"] = mean_hit_thr
report_df["mean_3_5sd_sum_threshold"] = mean_sum_thr
report_df["mean_3_5sd_sum_metric"] = mean_sum_log_hits_col
report_df["mean_3_5sd_threshold_source"] = mean_thresholds["threshold_source_label"]
report_df["q100_plus_median_hit_threshold"] = q100_hit_thr
report_df["q100_plus_median_sum_threshold"] = q100_sum_thr
report_df["q100_plus_median_sum_metric"] = q100med_sum_hits_col
report_df["q100_plus_median_threshold_source"] = q100_thresholds["threshold_source_label"]

# Binary z-score call based on thresholded mock-based z-scores
# default rule: reactive if >= 2 peptides have z > 3; indeterminate if exactly 1; otherwise non-reactive
def z_call_from_count(n):
    try:
        n = int(n)
    except Exception:
        return "Non-reactive"
    if n >= 2:
        return "Reactive"
    if n == 1:
        return "Indeterminate"
    return "Non-reactive"

report_df["call_z_mock"] = report_df["n_z_mock_above_3"].apply(z_call_from_count)
report_df["z_mock_binary_relaxed"] = report_df["call_z_mock"].apply(binary_call_relaxed)
report_df["z_mock_binary_strict"] = report_df["call_z_mock"].apply(binary_call_strict)

# Hybrid qualitative call: total signal intensity + peptide breadth.
if args.call_score == "mock_subtracted":
    # transform1 is the raw log10-CPM sum across the panel, so every panel peptide adds
    # roughly log10(background) to every sample alike. On a 50-peptide panel that shared
    # floor is ~60 points wide while a genuine response moves the score by a few points,
    # and t_high = negative median + 5 sigma can sit above the maximum achievable score.
    # The mock-subtracted transform (the same one --heatmap-mode ranked on) puts the
    # background near zero, and breadth is counted against the per-peptide mock threshold
    # rather than against the bare CPM floor, which every peptide clears in every well.
    hybrid_score_col = "transform_selected_method_continuous_sum"
    hybrid_hit_col = "n_above_mock_mean_3_5sd_panel"
else:
    hybrid_score_col = "transform1_raw_log10_cpm_sum"
    hybrid_hit_col = "panel_peptide_hit_count"

score_by_sample = dict(zip(report_df["sample"], report_df[hybrid_score_col]))
hit_by_sample = dict(zip(report_df["sample"], report_df[hybrid_hit_col]))
nonmock_for_thresholds = [s for s in sample_cols if s not in mock_cols]
hybrid_thresholds = derive_hybrid_call_thresholds(score_by_sample, hit_by_sample, neg_cols, nonmock_for_thresholds)
T_low = hybrid_thresholds["t_low"]
T_high = hybrid_thresholds["t_high"]

hybrid_min_negative_hits = int(args.hybrid_min_negative_hits) if args.hybrid_min_negative_hits is not None else 3
hybrid_min_positive_hits = int(args.hybrid_min_positive_hits) if args.hybrid_min_positive_hits is not None else int(derived_metric_min_hits)

report_df["call_hybrid"] = [
    hybrid_call_from_score_and_hits(
        score_by_sample.get(s, 0.0),
        hit_by_sample.get(s, 0),
        T_low,
        T_high,
        unsupervised_mode=unsupervised_mode,
        h_min_negative=hybrid_min_negative_hits,
        h_min_positive=hybrid_min_positive_hits
    )
    for s in report_df["sample"]
]
report_df["hybrid_score_threshold_low"] = T_low
report_df["hybrid_score_threshold_high"] = T_high
report_df["hybrid_threshold_source"] = hybrid_thresholds["source_label"]
report_df["hybrid_score_metric"] = hybrid_score_col
report_df["hybrid_hit_metric"] = hybrid_hit_col

call_mode_choice = (args.call_mode or (preset["call_mode"] if preset is not None else "hybrid")).strip().lower()
if call_mode_choice not in {"mean_3_5sd", "q100_plus_median", "hybrid", "both", "all"}:
    print(f"WARNING: Unrecognised call mode '{call_mode_choice}', defaulting to hybrid.")
    call_mode_choice = "hybrid"

primary_call_col = {
    "mean_3_5sd": "call_mean_3_5sd",
    "q100_plus_median": "call_q100_plus_median",
    "hybrid": "call_hybrid",
    "both": "call_hybrid",
    "all": "call_hybrid",
}[call_mode_choice]
report_df["call_primary"] = report_df[primary_call_col]
report_df["call"] = report_df["call_primary"]

report_df["call_mean_3_5sd_binary_relaxed"] = report_df["call_mean_3_5sd"].apply(binary_call_relaxed)
report_df["call_mean_3_5sd_binary_strict"] = report_df["call_mean_3_5sd"].apply(binary_call_strict)
report_df["call_q100_plus_median_binary_relaxed"] = report_df["call_q100_plus_median"].apply(binary_call_relaxed)
report_df["call_q100_plus_median_binary_strict"] = report_df["call_q100_plus_median"].apply(binary_call_strict)
report_df["call_hybrid_binary_relaxed"] = report_df["call_hybrid"].apply(binary_call_relaxed)
report_df["call_hybrid_binary_strict"] = report_df["call_hybrid"].apply(binary_call_strict)
report_df["call_binary_relaxed"] = report_df["call"].apply(binary_call_relaxed)
report_df["call_binary_strict"] = report_df["call"].apply(binary_call_strict)

report_df.to_csv(report_out, index=False)
print(f"\nmean_3_5sd call thresholds from {mean_thresholds['threshold_source_label']}:")
print(f"  hit threshold = {mean_hit_thr}")
print(f"  hits-only log10(CPM+1) sum threshold = {round(mean_sum_thr, 3)}")
print(f"\nq100_plus_median call thresholds from {q100_thresholds['threshold_source_label']}:")
print(f"  hit threshold = {q100_hit_thr}")
print(f"  hits-only sumCPM threshold = {round(q100_sum_thr, 3)}")
print(f"\nHybrid call thresholds from {hybrid_thresholds['source_label']}:")
print(f"  T_low  = {round(T_low, 3)}")
print(f"  T_high = {round(T_high, 3)}")
print(f"  Hit-count rule: Negative if <{hybrid_min_negative_hits} hits; Positive if >= {hybrid_min_positive_hits} hits and score >= T_high; otherwise Indeterminate when score >= T_low.")
if hybrid_thresholds.get("warning"):
    print(f"WARNING:{hybrid_thresholds['warning']}")
print(f"Primary call mode for generic call columns: {call_mode_choice} ({primary_call_col})")
print(f"\nWrote seropositivity / reactivity report: {report_out}")

elisa_merged_out = None
elisa_summary_out = None
elisa_stats_df = None
roc_plot_out = None
roc_auc_out = None
roc_points_out = None
orf_profile_plot_out = None
orf_profile_csv_out = None
if args.elisa_file:
    elisa_df = load_elisa_binary_csv(args.elisa_file)

    report_df = report_df.sort_values("sample").reset_index(drop=True).copy()
    report_df["sample_merge"] = report_df["sample"].apply(normalise_sample_for_merge)

    merged = (
        report_df.merge(
            elisa_df[["sample_merge", "elisa_binary"]],
            on="sample_merge",
            how="left"
        )
        .sort_values("sample")
        .reset_index(drop=True)
    )

    n_matched = int(merged["elisa_binary"].notna().sum())
    print(f"ELISA matches: {n_matched}/{len(merged)} samples")
    if n_matched == 0:
        print("WARNING: No ELISA values matched samples — check sample names and prefixes")

    # drop merge key from exported table
    merged_export = merged.drop(columns=["sample_merge"], errors="ignore")

    elisa_merged_out = f"{virus_prefix}_seropositivity_report_with_elisa.csv"
    merged_export.to_csv(elisa_merged_out, index=False)
    print(f"Wrote merged report with ELISA: {elisa_merged_out}")

    elisa_summary_out = f"{virus_prefix}_elisa_stats.csv"

    exclude_from_final_stats = get_samples_to_exclude_from_final_stats(selected_pos_cols, selected_neg_cols)

    elisa_stats_df = summarise_elisa_relationships(merged, elisa_summary_out, exclude_samples=exclude_from_final_stats)
    print(f"Wrote ELISA stats table: {elisa_summary_out}")
    roc_plot_out, roc_auc_out, roc_points_out = write_roc_outputs(merged, virus_prefix, exclude_samples=exclude_from_final_stats)
    orf_profile_res = write_orf_profile_plot(df, sample_cols, virus_prefix, merged)
    if orf_profile_res is not None:
        orf_profile_plot_out, orf_profile_csv_out = orf_profile_res
        print(f"Wrote ORF profile plot: {orf_profile_plot_out}")
        print(f"Wrote ORF profile summary table: {orf_profile_csv_out}")
    print(f"Wrote ROC plot: {roc_plot_out}")
    print(f"Wrote ROC AUC table: {roc_auc_out}")
    print(f"Wrote ROC points table: {roc_points_out}")
    print(f"Included continuous metrics: raw logCPM+1; mock median-corrected logCPM+1; other mock-corrected transforms; selected-method continuous metric ({selected_continuous_metric}); min{derived_metric_min_hits}-gated and hit-weighted raw logCPM+1 metrics")
    if "exclude_from_final_stats" in locals() and exclude_from_final_stats:
        print("Excluded selected comparison/negative controls from final ELISA stats:")
        for _s in sorted(exclude_from_final_stats):
            print(f"  {_s}")

    if elisa_stats_df is not None and not elisa_stats_df.empty:
        print("\n=== ELISA STATISTICS SUMMARY ===")
        for _, row in elisa_stats_df.iterrows():
            stat = row["statistic"]
            pval = row["p_value"]
            try:
                stat_txt = f"{float(stat):.4g}"
            except Exception:
                stat_txt = str(stat)
            try:
                p_txt = f"{float(pval):.4g}"
            except Exception:
                p_txt = str(pval)
            print(f"{row['comparison']}: statistic={stat_txt}, p={p_txt}")
        print(f"Full ELISA stats written to: {elisa_summary_out}\n")

panel_out = f"{virus_prefix}_selected_panel.txt"
with open(panel_out, "w") as f:
    for pep in final["panel_peptides"]:
        f.write(pep + "\n")
print(f"Wrote selected peptide panel: {panel_out}")

# ------------------------------------------------------------
# heatmap using selected panel
# ------------------------------------------------------------
panel_peptides = final["panel_peptides"]
heat_df = df[df["peptide"].isin(panel_peptides)].copy()

panel_order = {p: i for i, p in enumerate(panel_peptides)}
heat_df["panel_order"] = heat_df["peptide"].map(panel_order)
heat_df = heat_df.sort_values("panel_order").reset_index(drop=True)

ordered_samples = []
if pos_cols:
    ordered_samples.extend([c for c in pos_cols if c in sample_cols])
ordered_samples.extend([c for c in neg_cols if c in sample_cols and c not in ordered_samples])
ordered_samples.extend([c for c in sample_cols if c not in ordered_samples and c not in mock_cols])
ordered_samples.extend([c for c in mock_cols if c not in ordered_samples])

heat_raw = heat_df.set_index("peptide")[ordered_samples].copy().astype(float)
full_peptide_names = heat_raw.index.tolist()


call_mode_choice = (args.call_mode or (preset["call_mode"] if preset is not None else "hybrid")).strip().lower()

if args.heatmap_mode is not None:
    heatmap_mode = str(args.heatmap_mode).strip().lower()
elif preset is not None:
    heatmap_mode = str(preset["heatmap_mode"]).strip().lower()
else:
    default_heatmap_prompt = "8" if str(call_mode_choice).strip().lower() == "mean_3_5sd" else "7"
    heatmap_mode = prompt_with_default(
        "Choose heatmap mode: 0=no heatmap, 1=raw CPM (log10(CPM+1)), 2=z scores versus mocks, 3=log10(max(0, CPM - per-peptide mock mean - residual floor) + 1), 4=log10(max(0, CPM - per-peptide positive-mock median - residual floor) + 1), 5=log10(max(0, CPM - per-peptide mock Q75 - residual floor) + 1), 6=log10(max(0, CPM - per-peptide mock Q100/max - residual floor) + 1), 7=log10(max(0, CPM - per-peptide mock Q100/max - per-peptide mock median - residual floor) + 1), 8=log10(CPM+1) but only for values above mean + 3.5 SD threshold, 9=all heatmaps",
        default_heatmap_prompt
    ).strip().lower()

mode_aliases = {
    "0": "none",
    "none": "none",
    "no": "none",
    "no_heatmap": "none",
    "1": "raw_cpm",
    "raw": "raw_cpm",
    "raw_cpm": "raw_cpm",
    "2": "zscore",
    "z": "zscore",
    "zscore": "zscore",
    "z_scores": "zscore",
    "3": "minus_mean",
    "minus_mean": "minus_mean",
    "minus_mock_mean": "minus_mean",
    "4": "positive_mock_median",
    "positive_mock_median": "positive_mock_median",
    "5": "mock_q75",
    "q75": "mock_q75",
    "mock_q75": "mock_q75",
    "6": "mock_q100",
    "q100": "mock_q100",
    "max": "mock_q100",
    "mock_q100": "mock_q100",
    "7": "mock_q100_plus_median",
    "q100_plus_median": "mock_q100_plus_median",
    "mock_q100_plus_median": "mock_q100_plus_median",
    "8": "mean_3_5sd",
    "mean_3_5sd": "mean_3_5sd",
    "mean_3.5sd": "mean_3_5sd",
    "mean_heatmap": "mean_3_5sd",
    "9": "all",
    "all": "all",
}
heatmap_mode = mode_aliases.get(heatmap_mode, "mock_q100_plus_median")

# Modes to generate when --heatmap-mode all is selected.
# These suffixes must match the names used by build_heatmap(), because downstream
# scripts, including plot_peptide_seroprevalence.R in panel mode, look for files
# named <prefix>_heatmap_positive_values_<method>.csv.
all_modes = [
    "raw_cpm",
    "zscore",
    "minus_mean",
    "positive_mock_median",
    "mock_q75",
    "mock_q100",
    "mock_q100_plus_median",
    "mean_3_5sd",
]

if heatmap_mode == "none":
    include_mocks_in_heatmap = False
    if shared_positive_peptides is not None:
        min_positive_peptides_for_display = int(shared_positive_peptides)
    elif args.min_positive_peptides is not None:
        min_positive_peptides_for_display = int(args.min_positive_peptides)
    elif preset is not None:
        min_positive_peptides_for_display = int(preset["min_positive_peptides_for_display"])
    else:
        min_positive_peptides_for_display = 5
    z_clip = None
    displayed_samples = []
    heat_raw_display = pd.DataFrame()
    print("No heatmap will be produced.")
else:
    if (args.include_mocks or args.mocks) and args.exclude_mocks:
        parser.error("Use either --include-mocks/--mocks or --exclude-mocks, not both.")
    if args.include_mocks or args.mocks:
        include_mocks_in_heatmap = True
    elif args.exclude_mocks:
        include_mocks_in_heatmap = False
    elif preset is not None:
        include_mocks_in_heatmap = bool(preset["include_mocks_in_heatmap"])
    else:
        include_mocks_in_heatmap = yes_no_prompt(
            "Include mock columns in the heatmap display?",
            default="y"
        )

    if shared_positive_peptides is not None:
        min_positive_peptides_for_display = int(shared_positive_peptides)
    elif args.min_positive_peptides is not None:
        min_positive_peptides_for_display = int(args.min_positive_peptides)
    elif preset is not None:
        min_positive_peptides_for_display = int(preset["min_positive_peptides_for_display"])
    else:
        min_positive_peptides_for_display = int(prompt_with_default(
            "Minimum number of positive peptides required to display a sample column (0 disables suppression)",
            "5"
        ))

    z_clip = None
    if heatmap_mode in {"zscore", "all"}:
        if args.z_clip is not None:
            z_clip = float(args.z_clip)
        elif preset is not None:
            z_clip = float(preset["z_clip"])
        else:
            z_clip = float(prompt_with_default(
                "Enter the z-score display cutoff (values will be clipped to 0 to this value)",
                str(STRICT_SD_MULTIPLIER)
            ))

    if shared_positive_peptides is not None:
        print(f"Using shared positive-peptide threshold for calling and heatmap sample suppression: {shared_positive_peptides}")

    displayed_samples = ordered_samples[:] if include_mocks_in_heatmap else [c for c in ordered_samples if c not in mock_cols]

    if not displayed_samples:
        print("ERROR: No columns remain for heatmap display after excluding mocks.")
        sys.exit(1)

    heat_raw_display = heat_raw[displayed_samples].copy()

def apply_sample_suppression(heat_df_matrix, min_k):
    if min_k <= 0:
        return heat_df_matrix
    positive_counts = (heat_df_matrix.values > 0).sum(axis=0)
    keep_mask = positive_counts >= min_k
    suppressed = heat_df_matrix.copy()
    suppressed.loc[:, ~keep_mask] = 0.0
    return suppressed


def build_heatmap(mode_name):
    residual_floor = float(cpm_choice)
    floor_label = int(residual_floor) if residual_floor.is_integer() else residual_floor

    index = heat_raw_display.index.tolist()
    cols = heat_raw_display.columns.tolist()
    vals_out = np.zeros((len(index), len(cols)), dtype=float)
    heat_lookup = heat_df.set_index("peptide")

    for i, pep in enumerate(index):
        row = heat_lookup.loc[pep]
        for j, s in enumerate(cols):
            v = pd.to_numeric(row.get(s, np.nan), errors="coerce")
            vv = 0.0 if pd.isna(v) else float(v)
            stats = get_linked_mock_stats_from_row(row, s, mock_cols)

            if mode_name == "raw_cpm":
                vals_out[i, j] = math.log10(vv + 1.0)
            elif mode_name == "zscore":
                if stats["sd"] > 0:
                    z = (vv - stats["mean"]) / stats["sd"]
                else:
                    z = z_clip if vv > stats["mean"] else 0.0
                vals_out[i, j] = min(max(z, 0.0), z_clip)
            elif mode_name == "minus_mean":
                vals_out[i, j] = math.log10(max(0.0, vv - stats["mean"] - residual_floor) + 1.0)
            elif mode_name == "positive_mock_median":
                vals_out[i, j] = math.log10(max(0.0, vv - stats["positive_median"] - residual_floor) + 1.0)
            elif mode_name == "mock_q75":
                vals_out[i, j] = math.log10(max(0.0, vv - stats["q75"] - residual_floor) + 1.0)
            elif mode_name == "mock_q100":
                vals_out[i, j] = math.log10(max(0.0, vv - stats["q100"] - residual_floor) + 1.0)
            elif mode_name == "mock_q100_plus_median":
                vals_out[i, j] = math.log10(max(0.0, vv - stats["q100"] - stats["median"] - residual_floor) + 1.0)
            elif mode_name == "mean_3_5sd":
                thr = max(float(cpm_choice), stats["mean"] + STRICT_SD_MULTIPLIER * stats["sd"])
                vals_out[i, j] = math.log10(vv + 1.0) if vv > thr else 0.0
            else:
                raise ValueError(f"Unknown heatmap mode: {mode_name}")

    heat = pd.DataFrame(vals_out, index=index, columns=cols)
    heat = apply_sample_suppression(heat, min_positive_peptides_for_display)

    if mode_name == "raw_cpm":
        title_mode = "selected peptide panel heatmap (raw CPM)"
        cbar_label = "log10(CPM + 1)"
        heatmap_suffix = "raw_cpm"
    elif mode_name == "zscore":
        title_mode = f"selected peptide panel heatmap (z scores vs linked mocks; clipped at 0–{z_clip})"
        cbar_label = "z score relative to linked mocks (>=0 only)"
        heatmap_suffix = "zscore"
    elif mode_name == "minus_mean":
        title_mode = f"selected peptide panel heatmap (log10 of max(0, CPM - linked mock mean - residual floor {floor_label}))"
        cbar_label = f"log10(max(0, CPM - linked mock mean - {floor_label}) + 1)"
        heatmap_suffix = "minus_mean"
    elif mode_name == "positive_mock_median":
        title_mode = f"selected peptide panel heatmap (log10 of max(0, CPM - linked positive-mock median - residual floor {floor_label}))"
        cbar_label = f"log10(max(0, CPM - linked positive-mock median - {floor_label}) + 1)"
        heatmap_suffix = "positive_mock_median"
    elif mode_name == "mock_q75":
        title_mode = f"selected peptide panel heatmap (log10 of max(0, CPM - linked mock Q75 - residual floor {floor_label}))"
        cbar_label = f"log10(max(0, CPM - linked mock Q75 - {floor_label}) + 1)"
        heatmap_suffix = "mock_q75"
    elif mode_name == "mock_q100":
        title_mode = f"selected peptide panel heatmap (log10 of max(0, CPM - linked mock Q100/max - residual floor {floor_label}))"
        cbar_label = f"log10(max(0, CPM - linked mock Q100/max - {floor_label}) + 1)"
        heatmap_suffix = "mock_q100"
    elif mode_name == "mock_q100_plus_median":
        title_mode = f"selected peptide panel heatmap (log10 of max(0, CPM - linked mock Q100/max - linked mock median - residual floor {floor_label}))"
        cbar_label = f"log10(max(0, CPM - linked mock Q100/max - linked mock median - {floor_label}) + 1)"
        heatmap_suffix = "mock_q100_plus_median"
    elif mode_name == "mean_3_5sd":
        title_mode = "selected peptide panel heatmap (log10(CPM+1) only for values above linked mean + 3.5 SD)"
        cbar_label = "log10(CPM + 1)"
        heatmap_suffix = "mean_3_5sd"

    full_names_for_labels = [str(p) for p in heat.index.tolist()]
    plot_names = make_unique(full_names_for_labels)
    heat.index = plot_names

    fig, ax = plt.subplots(figsize=(max(12, len(heat.columns) * 0.35), max(6, len(heat.index) * 0.18)))
    im = ax.imshow(heat.values.astype(float), aspect='auto', interpolation='nearest')
    ax.set_xticks(range(len(heat.columns)))
    ax.set_xticklabels(heat.columns.tolist(), rotation=90, fontsize=8)
    ax.set_yticks(range(len(heat.index)))
    ax.set_yticklabels(heat.index.tolist(), fontsize=6)
    title_extra = ""
    if min_positive_peptides_for_display > 0:
        title_extra = f"; sample columns with <{min_positive_peptides_for_display} positive peptides suppressed"
    ax.set_title(f"{virus_prefix}: {title_mode}{title_extra}")
    ax.set_xlabel("Samples")
    ax.set_ylabel("Selected panel peptides")
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(cbar_label)
    plt.tight_layout()
    heatmap_out = f"{virus_prefix}_heatmap_{heatmap_suffix}.png"
    plt.savefig(heatmap_out, dpi=300)
    plt.close(fig)
    print(f"Wrote heatmap: {heatmap_out}")

    labels_out = None
    if args.output_set == "full":
        labels_out = f"{virus_prefix}_heatmap_labels_{heatmap_suffix}.csv"
        pd.DataFrame({"plot_label": plot_names, "full_peptide_name": full_names_for_labels}).to_csv(labels_out, index=False)
        print(f"Wrote heatmap label file: {labels_out}")

    # Write a full selected-panel signal matrix for downstream seroprevalence plotting.
    # Earlier versions dropped all-zero rows/columns, which made downstream group
    # denominators ambiguous unless raw_cpm was also generated. Keeping the full
    # displayed sample set allows the wrapper to generate only the chosen heatmap mode.
    positive_matrix = heat.copy()
    positives_out = f"{virus_prefix}_heatmap_positive_values_{heatmap_suffix}.csv"
    positive_matrix.to_csv(positives_out)
    print(f"Wrote full signal-value table: {positives_out}")

    return heatmap_out, labels_out, positives_out

heatmap_outputs = []
if heatmap_mode == "none":
    heatmap_outputs = []
elif heatmap_mode == "all":
    for mode_name in all_modes:
        heatmap_outputs.append(build_heatmap(mode_name))
else:
    heatmap_outputs.append(build_heatmap(heatmap_mode))

print("\nPreview of final calls:")
preview_df = report_df.sort_values(["group", "sample"])
for _, r in preview_df.iterrows():
    print(
        f"{r['sample']}	{r['group']}	"
        f"hits={int(r['n_above_mock_mean_3_5sd_panel'])}	"
        f"sumCPM={round(float(r['sum_cpm_panel']), 2)}	"
        f"{r['call']}"
    )

print("\nDone.")
print("Final outputs:")
if 'filtered_matrix_out' in locals():
    print(f"  {filtered_matrix_out}")
print(f"  {ranked_out}")
print(f"  {tuning_out}")
print(f"  {panel_out}")
print(f"  {report_out}")
if 'settings_out' in locals():
    print(f"  {settings_out}")
for heatmap_out, labels_out, positives_out in heatmap_outputs:
    print(f"  {heatmap_out}")
    if labels_out is not None:
        print(f"  {labels_out}")
    print(f"  {positives_out}")
