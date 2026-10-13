#!/usr/bin/env bash
# =============================================================================
# Launch the browser GUI for the PhIP-seq rank/call/plot pipeline v1.7.2.
#
# Same runs as run_phipseq_rank_call_and_plot_v1.7.2.sh, with every prompt
# answered on one page before the run starts, and the figures shown in the
# browser as they are written.
#
#   ./run_phipseq_gui.sh              # http://127.0.0.1:8765
#   ./run_phipseq_gui.sh 9000         # a different port
#
# Working on a remote machine? Forward the port from your laptop:
#   ssh -L 8765:127.0.0.1:8765 user@this-host
# then open http://127.0.0.1:8765 locally.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PORT="${1:-8765}"

exec python3 "${SCRIPT_DIR}/phipseq_gui.py" --port "${PORT}" "${@:2}"
