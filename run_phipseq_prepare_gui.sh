#!/usr/bin/env bash
# =============================================================================
# Launch the remote/SSH GUI for virus preparation + the rank/call pipeline.
#
#   ./run_phipseq_prepare_gui.sh              # http://127.0.0.1:8770
#   ./run_phipseq_prepare_gui.sh 9000         # a different port
#
# Sign in with your own Linux account on the analysis server. The password is
# used once to open the SSH connection and then discarded.
#
# Colleagues on other machines should tunnel rather than exposing the port:
#   ssh -L 8770:127.0.0.1:8770 <host-running-this>
# then open http://127.0.0.1:8770 locally. Credentials are refused from
# non-local addresses unless --allow-insecure-credentials is passed, because
# plain HTTP would put a real account password on the wire.
# =============================================================================
# ---------------------------------------------------------------- site config
# The only place this install's own paths appear. A colleague copying the
# folder to another server edits these three lines and nothing else; the Python
# is untouched. Each is a default for the login form, not a restriction, and
# each can be overridden in the environment before calling this script.
#
#   PHIPSEQ_SERVER    analysis server the GUI logs in to, prefilled on the
#                     login page.
#   PHIPSEQ_ROOT      project area a session is confined to. Must be wide
#                     enough to hold BOTH the input data and the results.
#   PHIPSEQ_PYTHON    interpreter on the analysis server. Needs pandas,
#                     numpy, matplotlib and scipy, and BLAST+ on its PATH.
#   PHIPSEQ_PIPELINE_DIR  left unset on purpose: it then defaults to the
#                     folder this script sits in, which is correct whenever
#                     the GUI travels with the scripts it drives.
# This install's own values live in site.env, which is NOT tracked by git, so
# nobody's hostname or directory layout is published and nobody has to undo
# someone else's defaults. Copy site.env.example to site.env and edit it.
SITE_ENV="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/site.env"
# shellcheck disable=SC1090
[ -f "$SITE_ENV" ] && . "$SITE_ENV"
export PHIPSEQ_SERVER="${PHIPSEQ_SERVER:-}"
export PHIPSEQ_ROOT="${PHIPSEQ_ROOT:-}"
export PHIPSEQ_PYTHON="${PHIPSEQ_PYTHON:-python3}"
export PHIPSEQ_REPO="${PHIPSEQ_REPO:-}"
# -----------------------------------------------------------------------------
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PORT="${1:-8770}"
exec python3 "${SCRIPT_DIR}/phipseq_prepare_gui.py" --port "${PORT}" "${@:2}"
