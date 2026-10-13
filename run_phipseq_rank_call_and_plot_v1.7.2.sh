#!/usr/bin/env bash
# =============================================================================
# PhIP-seq rank, call and plot — wrapper for pipeline v1.7.2
# -----------------------------------------------------------------------------
# Adapted from run_phipseq_rank_call_and_plot_v1.6.5_v3.sh, which hardcodes the
# v1.6.9 Python script. v1.6.x rejects --filter-term and --output-prefix, both
# introduced in v1.7.2, so that wrapper cannot drive these runs at all.
#
# Differences from the v1.6.5_v3 wrapper:
#   - targets phipseq_rank_call_and_plot_v1.7.2.py
#   - passes --filter-term and --output-prefix separately, so output filenames
#     no longer have to equal the peptide filter expression
#   - prompts for the output DIRECTORY and cd's there first. The pipeline writes
#     to the working directory, so running it from the wrong place silently
#     overwrites a previous run's outputs.
#   - prompts for --heatmap-mode, which in v1.7.2 also selects the peptide
#     RANKING method and therefore changes the panel and every downstream figure
#   - tees a run_stdout.log automatically, and runs python3 -u so that prompts
#     are not swallowed by pipe buffering
#   - drops the unused R_COVERAGE variable the old wrapper declared
#
# Mock and negative regexes are deliberately NOT prompted here: the Python
# script asks for them itself, and supplying them on the command line is what
# suppresses those questions.
# =============================================================================
set -euo pipefail

# A shell whose working directory has been deleted can still launch this script
# by absolute path, and the failure then surfaces much later as a Python
# FileNotFoundError on a bare output filename, which is a confusing way to
# learn that the cwd is gone. Fail here instead, with the actual fix.
if ! pwd -P >/dev/null 2>&1; then
  echo "ERROR: the current working directory no longer exists." >&2
  echo "       cd somewhere valid and re-run, e.g.:" >&2
  echo "       cd $(dirname "${BASH_SOURCE[0]}")" >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY_SCRIPT="${SCRIPT_DIR}/phipseq_rank_call_and_plot_v1.7.2.py"
R_COMPARE="${SCRIPT_DIR}/compare_runs_coverage.R"

if [[ ! -f "$PY_SCRIPT" ]]; then
  echo "ERROR: Cannot find Python script at: $PY_SCRIPT" >&2
  exit 1
fi

# ---------------------------------------------------------------- input matrix
if [[ -n "${PHIPSEQ_MATRIX:-}" && -f "${PHIPSEQ_MATRIX}" ]]; then
  MATRIX_CSV="${PHIPSEQ_MATRIX}"
  echo "Input CPM matrix CSV: ${MATRIX_CSV}"
else
  read -r -p "Path to input CPM matrix CSV: " MATRIX_CSV
  while [[ -z "${MATRIX_CSV}" || ! -f "${MATRIX_CSV}" ]]; do
    echo "Please enter a valid existing CSV file path."
    read -r -p "Path to input CPM matrix CSV: " MATRIX_CSV
  done
fi
MATRIX_CSV="$(cd "$(dirname "$MATRIX_CSV")" && pwd)/$(basename "$MATRIX_CSV")"

# ----------------------------------------------------------------- filter term
# Empty is not allowed: v1.7.2 treats a missing filter as "no filtering" and
# silently runs on the whole ~298k-peptide library rather than erroring.
read -r -p "Peptide filter term (e.g. 'Zaire_ebolavirus_GP' or 'Africa2 AND |s|'): " FILTER_TERM
while [[ -z "${FILTER_TERM}" ]]; do
  echo "A filter term is required — leaving it blank runs on the entire library."
  read -r -p "Peptide filter term: " FILTER_TERM
done

# --------------------------------------------------------------- output naming
read -r -p "Output prefix (used for every output filename): " OUTPUT_PREFIX
while [[ -z "${OUTPUT_PREFIX}" ]]; do
  read -r -p "Output prefix cannot be empty. Output prefix: " OUTPUT_PREFIX
done

read -r -p "Output directory (created if needed) [default: current directory]: " OUT_DIR
OUT_DIR="${OUT_DIR:-.}"
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"
echo "Outputs will be written to: $(pwd)"

if compgen -G "${OUTPUT_PREFIX}_*" > /dev/null; then
  read -r -p "WARNING: files matching '${OUTPUT_PREFIX}_*' already exist here and will be overwritten. Continue? [y/N]: " OVERWRITE
  [[ "${OVERWRITE,,}" =~ ^(y|yes)$ ]] || { echo "Aborted."; exit 1; }
fi

# --------------------------------------------------------------- heatmap mode
# In v1.7.2 --heatmap-mode also selects the peptide RANKING method, so it
# changes selected_panel.txt and every downstream figure, not just the heatmap.
#
# Modes 0 (no heatmap) and 9 (all heatmaps) are deliberately not offered: the
# ranking-side digit map (line ~1826) only knows 1-8, so 0 and 9 fall through
# to mock_q100_plus_median for RANKING while doing something else for the
# heatmap. Pass --heatmap-mode yourself at the pass-through prompt if you
# genuinely want one of those.
echo
echo "Heatmap / ranking mode — this sets BOTH the heatmap transform and the"
echo "peptide ranking method, so it determines which peptides enter the panel."
echo "  1  raw_cpm                log10(CPM+1), no background correction"
echo "  2  zscore                 z-scores versus mocks"
echo "  3  minus_mean             CPM minus per-peptide mock mean"
echo "  4  positive_mock_median   CPM minus per-peptide positive-mock median"
echo "  5  mock_q75               CPM minus per-peptide mock Q75      <- used by the exp23 runs"
echo "  6  mock_q100              CPM minus per-peptide mock Q100/max"
echo "  7  mock_q100_plus_median  Q100 + mock median                  <- pipeline default"
echo "  8  mean_3_5sd             log10(CPM+1) above mean + 3.5 SD only"
read -r -p "Mode (number or name) [default: 5 / mock_q75]: " HEATMAP_MODE
HEATMAP_MODE="${HEATMAP_MODE:-mock_q75}"

# Validate, because an unrecognised value is NOT an error in v1.7.2 — it
# silently becomes mock_q100_plus_median and quietly changes the panel.
while true; do
  case "${HEATMAP_MODE,,}" in
    1|raw_cpm|2|zscore|3|minus_mean|4|positive_mock_median|5|mock_q75|\
    6|mock_q100|7|mock_q100_plus_median|8|mean_3_5sd) break ;;
    *)
      echo "Unrecognised mode '${HEATMAP_MODE}' — it would silently fall back to mock_q100_plus_median."
      read -r -p "Mode (number or name) [default: 5 / mock_q75]: " HEATMAP_MODE
      HEATMAP_MODE="${HEATMAP_MODE:-mock_q75}"
      ;;
  esac
done

# --------------------------------------------------------------------- options
read -r -p "Path to ELISA CSV (optional; press Enter to skip): " ELISA_CSV
if [[ -n "${ELISA_CSV}" && ! -f "${ELISA_CSV}" ]]; then
  echo "WARNING: ELISA file not found; proceeding without it." >&2
  ELISA_CSV=""
fi

read -r -p "Use location collapse (identical locations only)? [y/N]: " USE_COLLAPSE
USE_COLLAPSE="${USE_COLLAPSE:-N}"

COLLAPSE_MODE=""
COLLAPSE_FLOOR=""
if [[ "${USE_COLLAPSE,,}" =~ ^(y|yes)$ ]]; then
  read -r -p "Collapse selection mode [raw_cpm/minus_mean_floor/minus_q75_floor/minus_q100_floor/minus_q100_plus_median_floor/zscore] (default raw_cpm): " COLLAPSE_MODE
  COLLAPSE_MODE="${COLLAPSE_MODE:-raw_cpm}"
  if [[ "${COLLAPSE_MODE}" != "raw_cpm" ]]; then
    read -r -p "Residual floor for mock-aware collapse (default 5): " COLLAPSE_FLOOR
    COLLAPSE_FLOOR="${COLLAPSE_FLOOR:-5}"
  fi
fi

# Word-split on whitespace, so no single pass-through argument may contain
# spaces. Anything needing spaces belongs in one of the prompts above.
#
# This is free text landing directly on the command line, and it follows a y/N
# question, so a stray "n" here is easy to type and comes back as an opaque
# "unrecognized arguments: n" from argparse. Require the first token to look
# like a flag, which every valid pass-through does.
echo
echo "Optional extra flags for the pipeline, e.g. '--panel-choice 20' or '--output-set full'."
read -r -p "Additional arguments (press Enter for none): " EXTRA_ARGS
while [[ -n "${EXTRA_ARGS// }" && "${EXTRA_ARGS}" != -* ]]; do
  echo "'${EXTRA_ARGS}' is not a flag — pipeline arguments start with '-'."
  echo "Press Enter to pass nothing."
  read -r -p "Additional arguments (press Enter for none): " EXTRA_ARGS
done

# ------------------------------------------------------------------------- run
CMD=(python3 -u "$PY_SCRIPT" "$MATRIX_CSV"
     --filter-term "$FILTER_TERM"
     --output-prefix "$OUTPUT_PREFIX"
     --heatmap-mode "$HEATMAP_MODE")

if [[ -n "${ELISA_CSV}" ]]; then
  CMD+=(--elisa-file "$ELISA_CSV")
fi
if [[ "${USE_COLLAPSE,,}" =~ ^(y|yes)$ ]]; then
  CMD+=(--collapse-by-location-max --collapse-selection-mode "$COLLAPSE_MODE")
  if [[ -n "$COLLAPSE_FLOOR" ]]; then
    CMD+=(--collapse-selection-floor "$COLLAPSE_FLOOR")
  fi
else
  CMD+=(--no-collapse-by-location-max)
fi
if [[ -n "${EXTRA_ARGS}" ]]; then
  # shellcheck disable=SC2206
  EXTRA_ARRAY=($EXTRA_ARGS)
  CMD+=("${EXTRA_ARRAY[@]}")
fi

echo
echo "Running:"
printf ' %q' "${CMD[@]}"
echo
echo

"${CMD[@]}" 2>&1 | tee run_stdout.log

# ------------------------------------------------- optional run-vs-run compare
echo
read -r -p "Compare with a second run directory? [y/N]: " RUN_COMPARE
RUN_COMPARE="${RUN_COMPARE:-N}"

if [[ "${RUN_COMPARE,,}" =~ ^(y|yes)$ ]]; then
  if [[ ! -f "$R_COMPARE" ]]; then
    echo "WARNING: Comparison script not found at: $R_COMPARE — skipping." >&2
  else
    read -r -p "Path to second run directory: " RUN2_DIR
    while [[ -z "${RUN2_DIR}" || ! -d "${RUN2_DIR}" ]]; do
      echo "Please enter a valid existing directory."
      read -r -p "Path to second run directory: " RUN2_DIR
    done
    DEFAULT_LABEL1="$(basename "$(pwd)")"
    DEFAULT_LABEL2="$(basename "${RUN2_DIR}")"
    read -r -p "Label for this run (default: ${DEFAULT_LABEL1}): " LABEL1
    LABEL1="${LABEL1:-${DEFAULT_LABEL1}}"
    read -r -p "Label for second run (default: ${DEFAULT_LABEL2}): " LABEL2
    LABEL2="${LABEL2:-${DEFAULT_LABEL2}}"
    echo "Running comparison plots..."
    Rscript "$R_COMPARE" "$(pwd)" "${RUN2_DIR}" "${LABEL1}" "${LABEL2}"
  fi
fi

echo
echo "All done. Log: $(pwd)/run_stdout.log"
