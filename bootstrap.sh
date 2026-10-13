#!/usr/bin/env bash
# =============================================================================
# Set up the PhIP-seq pipeline's dependencies on this machine, then verify them.
#
#   ./bootstrap.sh              # check; create the env only if something is missing
#   ./bootstrap.sh --check      # check and report, change nothing
#   ./bootstrap.sh --force      # recreate the environment even if it looks fine
#
# Safe to re-run. It never touches a working install, and it never needs root:
# everything lands inside your own conda installation.
#
# What it does, in order:
#   1. Looks for an interpreter that already satisfies the pipeline. If one
#      exists, it prints the path and stops -- nothing is downloaded.
#   2. Otherwise finds conda (or installs Miniforge into ~/miniforge3 if there
#      is none), and creates the 'phipseq' env from environment.yml.
#   3. Re-verifies, then prints the interpreter path to paste into the GUI's
#      "Python to run them with" field.
#
# WHY A SEPARATE STEP, AND NOT INSIDE THE ANALYSIS SCRIPTS
# -------------------------------------------------------
# If the interpreter lacks pandas, the fix is a *different interpreter*, and a
# running Python cannot become one. Dependencies therefore have to be settled
# before the pipeline starts, not during it.
#
# NOTE FOR CLUSTERS: run this on a login node. Compute nodes frequently have no
# outbound network, and the download will fail there in a confusing way. If your
# site already provides these as modules, you do not need this script at all --
# `module load` them and pass that interpreter to the GUI instead.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="${SCRIPT_DIR}/environment.yml"
ENV_NAME="phipseq"
MODULES=(pandas numpy matplotlib scipy)
BINARIES=(blastp tblastn makeblastdb)

MODE=create
case "${1:-}" in
  --check) MODE=check ;;
  --force) MODE=force ;;
  "") ;;
  *) echo "unknown option: $1 (expected --check, --force, or nothing)" >&2; exit 2 ;;
esac

say()  { printf '%s\n' "$*"; }
head2() { printf '\n== %s ==\n' "$*"; }

# --------------------------------------------------------------- verification
# Returns 0 only if this interpreter imports every module AND sees every binary.
# BLAST is checked through the interpreter's own PATH, because that is the PATH
# the pipeline's subprocess calls will inherit.
satisfies() {
  local py="$1"
  [ -x "$py" ] || command -v "$py" >/dev/null 2>&1 || return 1
  "$py" - <<'PYEOF' >/dev/null 2>&1 || return 1
import importlib, shutil, sys
for m in ("pandas", "numpy", "matplotlib", "scipy"):
    importlib.import_module(m)
missing = [b for b in ("blastp", "tblastn", "makeblastdb") if shutil.which(b) is None]
sys.exit(1 if missing else 0)
PYEOF
  return 0
}

report() {
  local py="$1"
  say "  interpreter : $py"
  "$py" - <<'PYEOF'
import importlib, shutil
for m in ("pandas", "numpy", "matplotlib", "scipy"):
    try:
        mod = importlib.import_module(m)
        print(f"  {m:<11s} {getattr(mod, '__version__', 'present')}")
    except ImportError:
        print(f"  {m:<11s} MISSING")
for b in ("blastp", "tblastn", "makeblastdb"):
    path = shutil.which(b)
    print(f"  {b:<12s} {path or 'MISSING'}")
PYEOF
}

# ------------------------------------------------- 1. is anything good already?
head2 "Looking for a usable interpreter"
CANDIDATES=()
[ -n "${PHIPSEQ_PYTHON:-}" ] && CANDIDATES+=("$PHIPSEQ_PYTHON")
if command -v conda >/dev/null 2>&1; then
  CONDA_BASE="$(conda info --base 2>/dev/null || true)"
  [ -n "$CONDA_BASE" ] && CANDIDATES+=("${CONDA_BASE}/envs/${ENV_NAME}/bin/python")
fi
CANDIDATES+=("$HOME/miniforge3/envs/${ENV_NAME}/bin/python"
             "$HOME/miniconda3/envs/${ENV_NAME}/bin/python"
             "$(command -v python3 || true)")

FOUND=""
SEEN=""
for py in "${CANDIDATES[@]}"; do
  [ -n "$py" ] || continue
  # The candidate list overlaps on purpose (conda's base may be ~/miniconda3),
  # so skip anything already tried rather than reporting it twice.
  case ":${SEEN}:" in *":${py}:"*) continue ;; esac
  SEEN="${SEEN}:${py}"
  if satisfies "$py"; then FOUND="$py"; break; fi
  say "  no  : $py"
done

if [ -n "$FOUND" ] && [ "$MODE" != force ]; then
  say "  yes : $FOUND"
  head2 "Already satisfied -- nothing downloaded"
  report "$FOUND"
  say ""
  say "Paste this into the GUI's \"Python to run them with\" field:"
  say "    $FOUND"
  exit 0
fi

if [ "$MODE" = check ]; then
  head2 "Not satisfied"
  say "No interpreter on this machine has all of: ${MODULES[*]} ${BINARIES[*]}"
  say "Re-run without --check to create the '${ENV_NAME}' environment."
  exit 1
fi

# ------------------------------------------------------------- 2. need conda?
head2 "Creating the '${ENV_NAME}' environment"
[ -f "$ENV_FILE" ] || { say "environment.yml not found beside this script"; exit 1; }

if ! command -v conda >/dev/null 2>&1; then
  say "No conda on this machine. Installing Miniforge into ~/miniforge3"
  say "(about 400 MB once the environment is built; no root required)."
  printf '  continue? [y/N] '
  read -r reply
  case "$reply" in
    [yY]*) ;;
    *) say "Stopped. Install conda yourself, then re-run this script."; exit 1 ;;
  esac
  TMP_INSTALLER="$(mktemp -t miniforge-XXXXXX.sh)"
  trap 'rm -f "$TMP_INSTALLER"' EXIT
  URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
  say "  downloading $URL"
  curl -fsSL "$URL" -o "$TMP_INSTALLER"
  bash "$TMP_INSTALLER" -b -p "$HOME/miniforge3"
  # shellcheck disable=SC1091
  . "$HOME/miniforge3/etc/profile.d/conda.sh"
fi

CONDA_BASE="$(conda info --base)"
# shellcheck disable=SC1091
. "${CONDA_BASE}/etc/profile.d/conda.sh"

if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  if [ "$MODE" = force ]; then
    say "  removing the existing '${ENV_NAME}' env (--force)"
    conda env remove -n "$ENV_NAME" -y
    conda env create -f "$ENV_FILE"
  else
    say "  '${ENV_NAME}' exists but is incomplete -- updating it"
    conda env update -n "$ENV_NAME" -f "$ENV_FILE" --prune
  fi
else
  conda env create -f "$ENV_FILE"
fi

# ------------------------------------------------------------- 3. verify again
ENV_PY="${CONDA_BASE}/envs/${ENV_NAME}/bin/python"
head2 "Verifying"
if satisfies "$ENV_PY"; then
  report "$ENV_PY"
  say ""
  say "Done. Paste this into the GUI's \"Python to run them with\" field:"
  say "    $ENV_PY"
  say ""
  say "For a shell instead:  conda activate ${ENV_NAME}"
else
  report "$ENV_PY"
  say ""
  say "The environment was created but still does not satisfy the pipeline."
  say "The lines marked MISSING above are what to chase."
  exit 1
fi
