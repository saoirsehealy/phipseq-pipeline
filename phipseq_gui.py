#!/usr/bin/env python3
# =============================================================================
# PhIP-seq rank, call and plot — browser GUI
# -----------------------------------------------------------------------------
# A local web front-end for run_phipseq_rank_call_and_plot_v1.7.2.sh.
#
# The shell wrapper asks its questions one at a time and the Python pipeline
# asks several more of its own part-way through the run, so a run cannot be
# reviewed before it starts and its figures only appear on disk afterwards.
# This GUI collects EVERY prompt from both layers on one page, shows the exact
# command it will run, then streams the log and displays the resulting plots
# and tables in the browser.
#
# How the prompts are satisfied
# -----------------------------
# Almost every interactive question in phipseq_rank_call_and_plot_v1.7.2.py has
# a matching command-line flag, so the GUI answers those by passing flags and
# the questions never appear. Two cases cannot be answered that way:
#
#   * "Do you have true negative controls on your run?" (line ~1625) has no
#     flag at all — --common-virus is parsed but never read by the script — so
#     the answer is fed on stdin. It is the first prompt the script can reach
#     once a collapse flag is supplied, which is why the GUI always passes one.
#   * Panel size and CPM cut-off may be left blank to accept the pipeline's own
#     tuned suggestion, which is only known mid-run. Blank fields simply omit
#     the flag; a blank line on stdin then accepts the printed default.
#
# stdin therefore carries: the true-negative-controls answer, then blank lines.
#
# Requires nothing outside the standard library. The pipeline itself is run as
# a subprocess with the interpreter given by --python (default python3), so it
# keeps using whatever environment already has pandas/matplotlib/scipy.
#
# Usage:
#   python3 phipseq_gui.py                 # http://127.0.0.1:8765
#   python3 phipseq_gui.py --port 9000 --open
# =============================================================================
import argparse
import csv
import json
import mimetypes
import os
import re
import shlex
import shutil
import subprocess
import sys
import threading
import time
import uuid
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse, parse_qs

HERE = os.path.dirname(os.path.abspath(__file__))
PY_SCRIPT = os.path.join(HERE, "phipseq_rank_call_and_plot_v1.7.2.py")
R_COMPARE = os.path.join(HERE, "compare_runs_coverage.R")
STATE_DIR = os.path.join(os.path.expanduser("~"), ".phipseq_gui")
LAST_FORM = os.path.join(STATE_DIR, "last_form.json")

PYTHON_BIN = "python3"

# ---------------------------------------------------------------------------
# Mirror of the pipeline's own sample-column detection, so the sample list the
# GUI shows carries the same 1-based indices the script would print for
# --pos-select / --neg-select. Kept in sync with lines 63-92 of the pipeline.
# ---------------------------------------------------------------------------
NON_SAMPLE_TRAILING_COLS = {
    "max_mock_cpm", "mean_mock_cpm", "prevalence_pct_nonmock",
    "ambiguous_score", "source_peptide_count", "orig_peptides",
}
DERIVED_SAMPLE_SUFFIXES = ("__ambiguous_score",)
DERIVED_SAMPLE_INFIXES = ("__mock__", "__max_mock_cpm", "__mean_mock_cpm", "__prevalence_pct_nonmock")


def is_true_sample_column(colname):
    col = str(colname)
    if col == "peptide" or col in NON_SAMPLE_TRAILING_COLS:
        return False
    if any(col.endswith(suf) for suf in DERIVED_SAMPLE_SUFFIXES):
        return False
    if any(infix in col for infix in DERIVED_SAMPLE_INFIXES):
        return False
    return True


# Modes 0 (no heatmap) and 9 (all heatmaps) are deliberately absent, for the
# reason the shell wrapper gives: the ranking-side digit map only knows 1-8, so
# 0 and 9 silently rank by mock_q100_plus_median while doing something else for
# the heatmap. Pass --heatmap-mode yourself under "additional arguments" if you
# genuinely want one of those.
HEATMAP_MODES = {
    "raw_cpm": "1 · raw_cpm — log10(CPM+1), no background correction",
    "zscore": "2 · zscore — z-scores versus mocks",
    "minus_mean": "3 · minus_mean — CPM minus per-peptide mock mean",
    "positive_mock_median": "4 · positive_mock_median — CPM minus per-peptide positive-mock median",
    "mock_q75": "5 · mock_q75 — CPM minus per-peptide mock Q75  (exp23 runs)",
    "mock_q100": "6 · mock_q100 — CPM minus per-peptide mock Q100/max",
    "mock_q100_plus_median": "7 · mock_q100_plus_median — Q100 + mock median  (pipeline default)",
    "mean_3_5sd": "8 · mean_3_5sd — log10(CPM+1) above mean + 3.5 SD only",
}
COLLAPSE_MODES = ["raw_cpm", "minus_mean_floor", "minus_q75_floor",
                  "minus_q100_floor", "minus_q100_plus_median_floor", "zscore"]
CALL_MODES = ["hybrid", "mean_3_5sd", "q100_plus_median", "both", "all"]

IMAGE_EXT = {".png", ".jpg", ".jpeg", ".svg", ".gif"}
TABLE_EXT = {".csv", ".tsv", ".txt", ".log"}


class FormError(Exception):
    """A problem with the submitted answers that the user must fix."""


def _s(form, key, default=""):
    value = form.get(key, default)
    return "" if value is None else str(value).strip()


def _abspath(path):
    return os.path.abspath(os.path.expanduser(path))


def build_command(form):
    """Turn the answered form into (cmd, stdin_answers, out_dir, notes).

    Raises FormError for anything the pipeline would only reject much later,
    or would quietly accept while doing something other than what was meant.
    """
    notes = []

    matrix = _s(form, "matrix_csv")
    if not matrix:
        raise FormError("Path to input CPM matrix CSV is required.")
    matrix = _abspath(matrix)
    if not os.path.isfile(matrix):
        raise FormError("Input CPM matrix CSV not found: %s" % matrix)

    # Empty is not allowed: v1.7.2 treats a missing filter as "no filtering"
    # and silently runs on the whole ~298k-peptide library rather than erroring.
    filter_term = _s(form, "filter_term")
    if not filter_term:
        raise FormError("A peptide filter term is required — leaving it blank runs on the entire library.")

    prefix = _s(form, "output_prefix")
    if not prefix:
        raise FormError("Output prefix is required; it names every output file.")
    if "/" in prefix or os.sep in prefix:
        raise FormError("Output prefix must be a filename fragment, not a path.")

    out_dir = _abspath(_s(form, "out_dir") or HERE)

    heatmap_mode = _s(form, "heatmap_mode", "mock_q75")
    if heatmap_mode not in HEATMAP_MODES:
        raise FormError("Unrecognised heatmap/ranking mode: %s" % heatmap_mode)

    cmd = [PYTHON_BIN, "-u", PY_SCRIPT, matrix,
           "--filter-term", filter_term,
           "--output-prefix", prefix,
           "--heatmap-mode", heatmap_mode]

    # ------------------------------------------------------------ collapse
    # Always pass one of the collapse flags. Besides being explicit, it keeps
    # the true-negative-controls question first among the prompts that remain,
    # which is what makes the stdin answer land where it is meant to.
    if _s(form, "collapse") in ("y", "yes", "on", "true", "1"):
        collapse_mode = _s(form, "collapse_mode", "raw_cpm") or "raw_cpm"
        if collapse_mode not in COLLAPSE_MODES:
            raise FormError("Unrecognised collapse selection mode: %s" % collapse_mode)
        cmd += ["--collapse-by-location-max", "--collapse-selection-mode", collapse_mode]
        if collapse_mode != "raw_cpm":
            floor = _s(form, "collapse_floor", "5") or "5"
            try:
                float(floor)
            except ValueError:
                raise FormError("Collapse residual floor must be a number: %s" % floor)
            cmd += ["--collapse-selection-floor", floor]
    else:
        cmd.append("--no-collapse-by-location-max")

    # -------------------------------------------------------- sample groups
    mock_regex = _s(form, "mock_regex", "mock") or "mock"
    cmd += ["--mock-regex", mock_regex]

    supervised = _s(form, "analysis_mode", "unsupervised") == "supervised"
    if supervised:
        pos_regex = _s(form, "pos_regex")
        neg_regex = _s(form, "neg_regex")
        pos_select = _s(form, "pos_select")
        neg_select = _s(form, "neg_select")
        if not pos_regex and not pos_select:
            raise FormError("Supervised mode needs known positives: give a positive regex "
                            "or positive sample indices.")
        has_true_negatives = _s(form, "true_negative_controls", "y") in ("y", "yes")
        if has_true_negatives and not neg_regex and not neg_select:
            raise FormError("Supervised mode needs known negatives unless you answer No to "
                            "\"true negative controls on this run\". Give a negative regex or "
                            "negative sample indices.")
        cmd.append("--supervised")
        if pos_regex:
            cmd += ["--pos-regex", pos_regex]
        if neg_regex:
            cmd += ["--neg-regex", neg_regex]
        if pos_select:
            cmd += ["--pos-select", pos_select]
        if neg_select:
            cmd += ["--neg-select", neg_select]
    else:
        # In unsupervised mode the script still asks for an optional
        # comparison regex used for labelling only. Passing --neg-regex, even
        # empty, is what suppresses that question.
        cmd += ["--unsupervised", "--neg-regex", _s(form, "comparison_regex")]

    # -------------------------------------------------------- panel tuning
    panel_sizes = _s(form, "panel_sizes", "5,10,15,20,25,30,40,50")
    if panel_sizes:
        for part in panel_sizes.split(","):
            if part.strip() and not part.strip().isdigit():
                raise FormError("Panel sizes must be whole numbers: %s" % panel_sizes)
        cmd += ["--panel-sizes", panel_sizes]

    cpm_thresholds = _s(form, "cpm_thresholds", "10,15,20,30")
    if cpm_thresholds:
        for part in cpm_thresholds.split(","):
            part = part.strip()
            if part:
                try:
                    float(part)
                except ValueError:
                    raise FormError("CPM cut-offs must be numbers: %s" % cpm_thresholds)
        cmd += ["--cpm-thresholds", cpm_thresholds]

    cpm_choice = _s(form, "cpm_choice")
    if cpm_choice:
        try:
            float(cpm_choice)
        except ValueError:
            raise FormError("CPM cut-off to use must be a number: %s" % cpm_choice)
        cmd += ["--cpm-choice", cpm_choice]
    else:
        notes.append("CPM cut-off left blank — the pipeline's own default "
                     "(the first CPM cut-off tested) is accepted on stdin.")

    panel_choice = _s(form, "panel_choice")
    if panel_choice:
        if not panel_choice.isdigit():
            raise FormError("Panel size to use must be a whole number: %s" % panel_choice)
        cmd += ["--panel-choice", panel_choice]
    else:
        notes.append("Panel size left blank — the pipeline's tuned suggestion is "
                     "accepted on stdin once tuning has run.")

    min_pep = _s(form, "min_positive_peptides", "5")
    if min_pep:
        if not min_pep.isdigit():
            raise FormError("Minimum positive peptides must be a whole number: %s" % min_pep)
        cmd += ["--min-positive-peptides", min_pep]

    # ---------------------------------------------------- calling / display
    call_mode = _s(form, "call_mode", "hybrid") or "hybrid"
    if call_mode not in CALL_MODES:
        raise FormError("Unrecognised call mode: %s" % call_mode)
    cmd += ["--call-mode", call_mode]

    if _s(form, "include_mocks", "y") in ("y", "yes", "on", "true", "1"):
        cmd.append("--include-mocks")
    else:
        cmd.append("--exclude-mocks")

    if heatmap_mode == "zscore":
        z_clip = _s(form, "z_clip", "3.5") or "3.5"
        try:
            float(z_clip)
        except ValueError:
            raise FormError("Z-score display cutoff must be a number: %s" % z_clip)
        cmd += ["--z-clip", z_clip]

    output_set = _s(form, "output_set", "minimal") or "minimal"
    if output_set not in ("minimal", "full"):
        raise FormError("Output set must be minimal or full.")
    cmd += ["--output-set", output_set]

    elisa = _s(form, "elisa_csv")
    if elisa:
        elisa = _abspath(elisa)
        if not os.path.isfile(elisa):
            raise FormError("ELISA CSV not found: %s" % elisa)
        cmd += ["--elisa-file", elisa]

    # ------------------------------------------------------- pass-through
    extra = _s(form, "extra_args")
    if extra:
        try:
            extra_tokens = shlex.split(extra)
        except ValueError as exc:
            raise FormError("Could not parse additional arguments: %s" % exc)
        if extra_tokens and not extra_tokens[0].startswith("-"):
            raise FormError("'%s' is not a flag — pipeline arguments start with '-'." % extra_tokens[0])
        cmd += extra_tokens

    # -------------------------------------------------------------- stdin
    # First answer is the true-negative-controls question, the only prompt with
    # no flag behind it. Everything after is a blank line, which every remaining
    # prompt reads as "use the printed default".
    tn = "n" if _s(form, "true_negative_controls", "y") in ("n", "no") else "y"
    stdin_answers = [tn] + [""] * 200
    if tn == "n":
        notes.append("No true negative controls: the run uses common-virus mode, "
                     "deriving final thresholds from mocks.")

    return cmd, stdin_answers, out_dir, notes


# ===========================================================================
# Run management
# ===========================================================================
class Run:
    def __init__(self, run_id, cmd, stdin_answers, out_dir, form, notes):
        self.id = run_id
        self.cmd = cmd
        self.stdin_answers = stdin_answers
        self.out_dir = out_dir
        # Directories this run is allowed to serve files from. compare_runs_coverage.R
        # writes its figures to dirname(run1_dir), not into the run directory, so
        # that parent is added as a second root when the comparison actually runs.
        self.roots = [out_dir]
        self.form = form
        self.notes = notes
        self.started = time.time()
        self.finished = None
        self.state = "starting"
        self.returncode = None
        self.error = None
        self.compare_note = None
        self.proc = None
        self.log = bytearray()
        self.lock = threading.Lock()

    # -- log ---------------------------------------------------------------
    def write(self, chunk):
        with self.lock:
            self.log.extend(chunk)

    def say(self, text):
        self.write(text.encode("utf-8", "replace"))

    def log_since(self, offset):
        with self.lock:
            total = len(self.log)
            if offset < 0 or offset > total:
                offset = 0
            return bytes(self.log[offset:total]).decode("utf-8", "replace"), total

    # -- outputs -----------------------------------------------------------
    def outputs(self):
        """Files this run created or refreshed, across all of its roots."""
        cutoff = self.started - 2
        images, tables, others = [], [], []
        for index, root in enumerate(self.roots):
            try:
                entries = list(os.scandir(root))
            except OSError:
                continue
            for entry in entries:
                if not entry.is_file():
                    continue
                # Beyond the run's own directory, claim only the comparison
                # figures — the parent holds other people's runs too.
                if index > 0 and not entry.name.startswith("comparison_"):
                    continue
                try:
                    stat = entry.stat()
                except OSError:
                    continue
                if stat.st_mtime < cutoff:
                    continue
                ext = os.path.splitext(entry.name)[1].lower()
                item = {"name": entry.name, "size": stat.st_size,
                        "mtime": stat.st_mtime, "root": index}
                if ext in IMAGE_EXT:
                    images.append(item)
                elif ext in TABLE_EXT:
                    tables.append(item)
                else:
                    others.append(item)
        images.sort(key=lambda i: i["name"])
        tables.sort(key=lambda i: i["name"])
        others.sort(key=lambda i: i["name"])
        return {"images": images, "tables": tables, "others": others}

    def status(self, offset=0):
        chunk, total = self.log_since(offset)
        return {
            "id": self.id,
            "state": self.state,
            "returncode": self.returncode,
            "error": self.error,
            "compare_note": self.compare_note,
            "out_dir": self.out_dir,
            "started": self.started,
            "finished": self.finished,
            "elapsed": (self.finished or time.time()) - self.started,
            "log": chunk,
            "offset": total,
            "outputs": self.outputs(),
            "notes": self.notes,
            "cmd": self.cmd,
        }

    # -- execution ---------------------------------------------------------
    def start(self):
        threading.Thread(target=self._run_guarded, daemon=True).start()

    def _pump(self, stream, log_fh):
        """Copy a child's output to the log as it arrives."""
        fd = stream.fileno()
        while True:
            try:
                chunk = os.read(fd, 65536)
            except OSError:
                break
            if not chunk:
                break
            self.write(chunk)
            if log_fh:
                try:
                    log_fh.write(chunk)
                    log_fh.flush()
                except OSError:
                    pass

    def _run_guarded(self):
        # Without this the run would sit at "running" for ever if the worker
        # thread died, with the browser polling a status that never changes.
        try:
            self._execute()
        except Exception as exc:
            self.state = "failed"
            self.error = "GUI error while running the pipeline: %s: %s" % (type(exc).__name__, exc)
            self.finished = time.time()
            self.say("\nGUI ERROR: %s: %s\n" % (type(exc).__name__, exc))
            if self.proc and self.proc.poll() is None:
                self.proc.kill()

    def _execute(self):
        log_path = os.path.join(self.out_dir, "run_stdout.log")
        try:
            os.makedirs(self.out_dir, exist_ok=True)
        except OSError as exc:
            self.state, self.error = "failed", "Could not create output directory: %s" % exc
            self.finished = time.time()
            return

        self.say("Output directory: %s\n\nRunning:\n  %s\n\n"
                 % (self.out_dir, " ".join(shlex.quote(c) for c in self.cmd)))
        self.say("Answers fed on stdin: %s, then blank lines (accept printed defaults).\n\n"
                 % self.stdin_answers[0])

        try:
            log_fh = open(log_path, "wb")
        except OSError as exc:
            log_fh = None
            self.say("WARNING: could not open %s for writing: %s\n" % (log_path, exc))

        try:
            self.proc = subprocess.Popen(
                self.cmd, cwd=self.out_dir,
                stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                bufsize=0,
            )
        except OSError as exc:
            self.state, self.error = "failed", "Could not start the pipeline: %s" % exc
            self.finished = time.time()
            if log_fh:
                log_fh.close()
            return

        self.state = "running"

        # Feed the answers up front. A couple of hundred short lines sit well
        # inside the pipe buffer, so this cannot deadlock against a child that
        # has not reached its first prompt yet.
        try:
            self.proc.stdin.write(("\n".join(self.stdin_answers) + "\n").encode())
            self.proc.stdin.flush()
            self.proc.stdin.close()
        except OSError:
            pass

        # Read raw bytes rather than lines: the pipeline writes its prompts
        # without a trailing newline, and line-buffered reads would hold them
        # back until the answer had already gone in. bufsize=0 means stdout is
        # an unbuffered raw stream, so read from the descriptor directly.
        self._pump(self.proc.stdout, log_fh)
        self.returncode = self.proc.wait()

        if self.returncode == 0:
            self._maybe_compare(log_fh)
            self.state = "finished"
        elif self.state == "stopping":
            self.state = "stopped"
        else:
            self.state = "failed"
            self.error = "The pipeline exited with code %s." % self.returncode

        self.finished = time.time()
        if log_fh:
            log_fh.close()

    def _maybe_compare(self, log_fh):
        form = self.form
        if _s(form, "compare_enabled") not in ("y", "yes", "on", "true", "1"):
            return
        run2 = _s(form, "compare_dir")
        if not run2 or not os.path.isdir(_abspath(run2)):
            self.compare_note = "Comparison skipped: second run directory not found."
            self.say("\n%s\n" % self.compare_note)
            return
        if not os.path.isfile(R_COMPARE):
            self.compare_note = "Comparison skipped: %s not found." % R_COMPARE
            self.say("\n%s\n" % self.compare_note)
            return
        if not shutil.which("Rscript"):
            self.compare_note = "Comparison skipped: Rscript is not on PATH."
            self.say("\n%s\n" % self.compare_note)
            return
        run2 = _abspath(run2)
        label1 = _s(form, "label1") or os.path.basename(self.out_dir)
        label2 = _s(form, "label2") or os.path.basename(run2)
        cmd = ["Rscript", R_COMPARE, self.out_dir, run2, label1, label2]
        self.say("\n" + "=" * 70 + "\nRunning comparison plots:\n  %s\n\n"
                 % " ".join(shlex.quote(c) for c in cmd))
        try:
            proc = subprocess.Popen(cmd, cwd=self.out_dir,
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=0)
        except OSError as exc:
            self.compare_note = "Comparison failed to start: %s" % exc
            self.say("%s\n" % self.compare_note)
            return
        self._pump(proc.stdout, log_fh)
        rc = proc.wait()
        if rc == 0:
            parent = os.path.dirname(self.out_dir)
            if parent and parent not in self.roots:
                self.roots.append(parent)
        self.compare_note = ("Comparison plots written to %s" % os.path.dirname(self.out_dir) if rc == 0
                             else "Comparison script exited with code %s — see the log." % rc)
        self.say("\n%s\n" % self.compare_note)

    def stop(self):
        if self.proc and self.proc.poll() is None:
            self.state = "stopping"
            self.proc.terminate()
            return True
        return False


RUNS = {}
RUNS_LOCK = threading.Lock()


# ===========================================================================
# Small filesystem helpers used by the browser UI
# ===========================================================================
def read_sample_columns(matrix_path, mock_regex):
    """Header-only read, so a 300k-row matrix still answers instantly."""
    path = _abspath(matrix_path)
    if not os.path.isfile(path):
        raise FormError("Matrix CSV not found: %s" % path)
    with open(path, "r", newline="") as fh:
        header = next(csv.reader(fh), [])
    if "peptide" not in header:
        raise FormError("The input matrix must contain a 'peptide' column.")
    samples = [c for c in header if is_true_sample_column(c)]
    try:
        pattern = re.compile(mock_regex or "mock")
    except re.error as exc:
        raise FormError("Mock regex is not a valid regular expression: %s" % exc)
    mocks = [c for c in samples if pattern.search(c)]
    nonmocks = [c for c in samples if c not in mocks]
    return {
        "n_columns": len(header),
        "n_samples": len(samples),
        "mocks": mocks,
        # 1-based indices, matching what --pos-select / --neg-select expect.
        "nonmocks": [{"index": i, "name": c} for i, c in enumerate(nonmocks, start=1)],
    }


def browse_dir(path, want):
    path = _abspath(path or HERE)
    if not os.path.isdir(path):
        path = HERE
    dirs, files = [], []
    try:
        for entry in os.scandir(path):
            if entry.name.startswith("."):
                continue
            if entry.is_dir():
                dirs.append(entry.name)
            elif want != "dir":
                ext = os.path.splitext(entry.name)[1].lower()
                if want == "csv" and ext not in (".csv", ".tsv", ".txt"):
                    continue
                files.append(entry.name)
    except OSError as exc:
        raise FormError("Cannot read directory: %s" % exc)
    dirs.sort(key=str.lower)
    files.sort(key=str.lower)
    parent = os.path.dirname(path.rstrip(os.sep)) or "/"
    return {"dir": path, "parent": parent, "dirs": dirs, "files": files}


def preview_table(path, limit=40):
    rows = []
    truncated = False
    ext = os.path.splitext(path)[1].lower()
    if ext in (".txt", ".log"):
        with open(path, "r", errors="replace") as fh:
            for i, line in enumerate(fh):
                if i >= limit:
                    truncated = True
                    break
                rows.append([line.rstrip("\n")])
        return {"header": None, "rows": rows, "truncated": truncated}
    with open(path, "r", newline="", errors="replace") as fh:
        reader = csv.reader(fh)
        header = next(reader, [])
        for i, row in enumerate(reader):
            if i >= limit:
                truncated = True
                break
            rows.append(row)
    return {"header": header, "rows": rows, "truncated": truncated}


def existing_prefix_files(out_dir, prefix):
    out_dir = _abspath(out_dir)
    if not os.path.isdir(out_dir):
        return []
    return sorted(n for n in os.listdir(out_dir) if n.startswith(prefix + "_"))


def save_last_form(form):
    try:
        os.makedirs(STATE_DIR, exist_ok=True)
        with open(LAST_FORM, "w") as fh:
            json.dump(form, fh, indent=2)
    except OSError:
        pass


def load_last_form():
    try:
        with open(LAST_FORM) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return {}


# ===========================================================================
# The page
# ===========================================================================
PAGE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>PhIP-seq rank, call and plot</title>
<style>
:root{
  --bg:#f6f7f9; --panel:#ffffff; --ink:#161a1f; --muted:#6b7480; --line:#dfe3e8;
  --accent:#2f6f4e; --accent-ink:#ffffff; --warn:#8a5a00; --warn-bg:#fff6e0;
  --err:#a4231f; --err-bg:#fdecec; --ok:#1f6b3a; --code:#f2f3f5;
}
@media (prefers-color-scheme: dark){
  :root{
    --bg:#14171b; --panel:#1b1f24; --ink:#e6e9ed; --muted:#9aa3ad; --line:#2c3238;
    --accent:#4e9c73; --accent-ink:#0d1116; --warn:#e0b054; --warn-bg:#332a15;
    --err:#f08a86; --err-bg:#3a1f1e; --ok:#7fd0a0; --code:#22272d;
  }
}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
  font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
header{position:sticky;top:0;z-index:20;display:flex;align-items:center;gap:14px;
  padding:10px 18px;background:var(--panel);border-bottom:1px solid var(--line)}
header h1{margin:0;font-size:15px;font-weight:650;letter-spacing:.01em}
header .sub{color:var(--muted);font-size:12px}
.pill{margin-left:auto;font-size:12px;padding:3px 10px;border-radius:999px;
  border:1px solid var(--line);color:var(--muted)}
.pill.running{border-color:var(--accent);color:var(--accent)}
.pill.finished{border-color:var(--ok);color:var(--ok)}
.pill.failed,.pill.stopped{border-color:var(--err);color:var(--err)}
main{display:grid;grid-template-columns:minmax(380px,1fr) minmax(420px,1.15fr);gap:14px;padding:14px;
  align-items:start}
@media (max-width:1050px){main{grid-template-columns:1fr}}
.panel{background:var(--panel);border:1px solid var(--line);border-radius:10px}
.col{display:flex;flex-direction:column;gap:14px;min-width:0}
fieldset{border:0;border-top:1px solid var(--line);margin:0;padding:14px 16px}
fieldset:first-of-type{border-top:0}
legend{padding:0;font-weight:650;font-size:12px;text-transform:uppercase;
  letter-spacing:.06em;color:var(--muted)}
.field{margin-top:12px}
.field > label{display:block;font-weight:600;margin-bottom:4px}
.hint{color:var(--muted);font-size:12px;margin-top:3px}
input[type=text],input[type=number],select,textarea{
  width:100%;padding:7px 9px;border:1px solid var(--line);border-radius:7px;
  background:var(--bg);color:var(--ink);font:inherit}
textarea{resize:vertical;min-height:58px}
input:focus,select:focus,textarea:focus{outline:2px solid var(--accent);outline-offset:-1px}
.row{display:flex;gap:8px;align-items:center}
.row > input{flex:1;min-width:0}
.grid2{display:grid;grid-template-columns:1fr 1fr;gap:12px}
button{font:inherit;padding:7px 13px;border-radius:7px;border:1px solid var(--line);
  background:var(--panel);color:var(--ink);cursor:pointer}
button:hover{border-color:var(--accent)}
button.primary{background:var(--accent);border-color:var(--accent);color:var(--accent-ink);font-weight:600}
button.small{padding:4px 9px;font-size:12px}
button:disabled{opacity:.5;cursor:not-allowed}
.actions{position:sticky;bottom:0;display:flex;gap:8px;align-items:center;flex-wrap:wrap;
  padding:12px 16px;background:var(--panel);border-top:1px solid var(--line);
  border-radius:0 0 10px 10px}
.msg{padding:10px 12px;border-radius:8px;margin:0 16px 14px;font-size:13px;white-space:pre-wrap}
.msg.err{background:var(--err-bg);color:var(--err)}
.msg.warn{background:var(--warn-bg);color:var(--warn)}
pre.cmd{margin:0 16px 14px;padding:10px 12px;background:var(--code);border-radius:8px;
  font-size:12px;white-space:pre-wrap;word-break:break-all}
.tabs{display:flex;gap:2px;padding:8px 10px 0;border-bottom:1px solid var(--line);flex-wrap:wrap}
.tabs button{border:1px solid transparent;border-bottom:0;border-radius:8px 8px 0 0;
  background:transparent;color:var(--muted);padding:6px 12px}
.tabs button.on{background:var(--bg);border-color:var(--line);color:var(--ink);font-weight:600}
.tabbody{padding:14px 16px;min-height:220px}
#log{margin:0;padding:12px;background:var(--code);border-radius:8px;height:60vh;overflow:auto;
  font:12px/1.45 ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;white-space:pre-wrap;word-break:break-word}
.gallery{display:flex;flex-direction:column;gap:16px}
.plot{border:1px solid var(--line);border-radius:9px;overflow:hidden;background:var(--bg)}
.plot .cap{display:flex;gap:10px;align-items:center;padding:7px 10px;border-bottom:1px solid var(--line);
  font-size:12px}
.plot .cap .name{font-weight:600;word-break:break-all}
.plot .cap .sz{color:var(--muted);margin-left:auto;white-space:nowrap}
.plot img{display:block;width:100%;height:auto;cursor:zoom-in;background:#fff}
.files{width:100%;border-collapse:collapse;font-size:13px}
.files td{padding:6px 8px;border-bottom:1px solid var(--line);vertical-align:top}
.files td.sz{color:var(--muted);text-align:right;white-space:nowrap}
.files a{color:var(--accent)}
.files .fname{cursor:pointer;color:var(--accent);word-break:break-all}
.tablewrap{overflow:auto;max-height:60vh;border:1px solid var(--line);border-radius:8px;margin-top:10px}
.tablewrap table{border-collapse:collapse;font-size:12px;white-space:nowrap}
.tablewrap th,.tablewrap td{padding:4px 9px;border-bottom:1px solid var(--line);text-align:left}
.tablewrap th{position:sticky;top:0;background:var(--panel);font-weight:650}
.empty{color:var(--muted);padding:24px 0;text-align:center}
.samplebox{margin-top:10px;border:1px solid var(--line);border-radius:8px;max-height:230px;overflow:auto}
.samplebox table{width:100%;border-collapse:collapse;font-size:12px}
.samplebox td{padding:3px 8px;border-bottom:1px solid var(--line)}
.samplebox td.i{color:var(--muted);width:38px;text-align:right}
.samplebox td.b{width:112px;white-space:nowrap}
.samplebox .mockrow{color:var(--muted)}
.modal{position:fixed;inset:0;background:rgba(0,0,0,.55);display:none;align-items:center;
  justify-content:center;z-index:50;padding:24px}
.modal.on{display:flex}
.modal .box{background:var(--panel);border-radius:10px;max-width:720px;width:100%;
  max-height:80vh;display:flex;flex-direction:column}
.modal .box header{position:static;border-radius:10px 10px 0 0}
.modal .list{overflow:auto;padding:6px 0}
.modal .list div{padding:6px 16px;cursor:pointer}
.modal .list div:hover{background:var(--bg)}
.modal .foot{padding:10px 16px;border-top:1px solid var(--line);display:flex;gap:8px}
#lightbox{position:fixed;inset:0;background:rgba(0,0,0,.9);display:none;z-index:60;
  overflow:auto;padding:20px;text-align:center}
#lightbox.on{display:block}
#lightbox img{max-width:none;background:#fff}
.tag{display:inline-block;font-size:11px;padding:1px 7px;border-radius:999px;
  border:1px solid var(--line);color:var(--muted);margin-left:6px}
</style>
</head>
<body>

<header>
  <h1>PhIP-seq rank, call and plot</h1>
  <span class="sub">pipeline v1.7.2</span>
  <span id="pill" class="pill">idle</span>
</header>

<main>
  <!-- ------------------------------------------------------ form column -->
  <div class="col">
  <form id="form" class="panel" onsubmit="return false">

    <fieldset>
      <legend>Input and output</legend>

      <div class="field">
        <label for="matrix_csv">Path to input CPM matrix CSV</label>
        <div class="row">
          <input type="text" id="matrix_csv" name="matrix_csv" placeholder="/path/to/matrix.csv">
          <button class="small" type="button" onclick="openBrowser('matrix_csv','csv')">Browse</button>
        </div>
      </div>

      <div class="field">
        <label for="filter_term">Peptide filter term</label>
        <input type="text" id="filter_term" name="filter_term"
               placeholder="Zaire_ebolavirus_GP    or    Africa2 AND |s|">
        <div class="hint">Required. Supports AND plus OR groups, e.g. <code>cchfv AND (|s| OR |m|) AND Africa3</code>.
          Blank is not allowed: v1.7.2 reads a missing filter as no filtering and runs the whole ~298k library.</div>
      </div>

      <div class="grid2">
        <div class="field">
          <label for="output_prefix">Output prefix</label>
          <input type="text" id="output_prefix" name="output_prefix" placeholder="zaire_GP_unsup">
          <div class="hint">Names every output file.</div>
        </div>
        <div class="field">
          <label for="out_dir">Output directory</label>
          <div class="row">
            <input type="text" id="out_dir" name="out_dir">
            <button class="small" type="button" onclick="openBrowser('out_dir','dir')">Browse</button>
          </div>
          <div class="hint">Created if needed. The pipeline writes to this directory.</div>
        </div>
      </div>
    </fieldset>

    <fieldset>
      <legend>Ranking and heatmap mode</legend>
      <div class="field">
        <label for="heatmap_mode">Mode</label>
        <select id="heatmap_mode" name="heatmap_mode"></select>
        <div class="hint">In v1.7.2 this sets <em>both</em> the heatmap transform and the peptide
          ranking method, so it decides which peptides enter the panel and changes every downstream
          figure. Modes 0 and 9 are not offered — their ranking silently falls back to
          mock_q100_plus_median.</div>
      </div>
    </fieldset>

    <fieldset>
      <legend>Sample groups</legend>

      <div class="grid2">
        <div class="field">
          <label for="mock_regex">Mock sample regex</label>
          <input type="text" id="mock_regex" name="mock_regex" value="mock">
        </div>
        <div class="field">
          <label for="true_negative_controls">True negative controls on this run?</label>
          <select id="true_negative_controls" name="true_negative_controls">
            <option value="y">Yes</option>
            <option value="n">No — common-virus mode (thresholds from mocks)</option>
          </select>
        </div>
      </div>

      <div class="field">
        <button class="small" type="button" onclick="loadColumns()">List sample columns from the matrix</button>
        <span id="colnote" class="hint" style="display:inline-block;margin-left:8px"></span>
        <div id="samplebox"></div>
      </div>

      <div class="field">
        <label for="analysis_mode">Analysis mode</label>
        <select id="analysis_mode" name="analysis_mode" onchange="syncMode()">
          <option value="unsupervised">Unsupervised — rank by non-mock enrichment above mocks</option>
          <option value="supervised">Supervised — use known positive and negative samples</option>
        </select>
      </div>

      <div id="unsup_fields">
        <div class="field">
          <label for="comparison_regex">Comparison / likely-negative regex <span class="tag">labelling only</span></label>
          <input type="text" id="comparison_regex" name="comparison_regex" placeholder="optional">
        </div>
      </div>

      <div id="sup_fields" style="display:none">
        <div class="grid2">
          <div class="field">
            <label for="pos_regex">Positive sample regex</label>
            <input type="text" id="pos_regex" name="pos_regex" placeholder="MHRA    or    AV|MHRA">
          </div>
          <div class="field">
            <label for="neg_regex">Negative sample regex</label>
            <input type="text" id="neg_regex" name="neg_regex" placeholder="HC17|KOVAG|PRIVIGEN">
          </div>
          <div class="field">
            <label for="pos_select">Positive sample indices</label>
            <input type="text" id="pos_select" name="pos_select" placeholder="4-6,9">
          </div>
          <div class="field">
            <label for="neg_select">Negative sample indices</label>
            <input type="text" id="neg_select" name="neg_select" placeholder="1-3,8">
          </div>
        </div>
        <div class="hint">Indices are 1-based over the non-mock samples listed above. Use the
          <strong>+pos</strong> / <strong>+neg</strong> buttons in that list to fill these in.</div>
      </div>
    </fieldset>

    <fieldset>
      <legend>Location collapse</legend>
      <div class="field">
        <label for="collapse">If several peptides share a location, keep only the maximum per sample?</label>
        <select id="collapse" name="collapse" onchange="syncCollapse()">
          <option value="n">No — keep all matched peptides</option>
          <option value="y">Yes — collapse to unique locations</option>
        </select>
      </div>
      <div id="collapse_fields" style="display:none" class="grid2">
        <div class="field">
          <label for="collapse_mode">Collapse selection mode</label>
          <select id="collapse_mode" name="collapse_mode" onchange="syncCollapse()"></select>
        </div>
        <div class="field" id="floor_field">
          <label for="collapse_floor">Residual floor</label>
          <input type="text" id="collapse_floor" name="collapse_floor" value="5">
        </div>
      </div>
    </fieldset>

    <fieldset>
      <legend>Panel tuning</legend>
      <div class="grid2">
        <div class="field">
          <label for="panel_sizes">Panel sizes to test</label>
          <input type="text" id="panel_sizes" name="panel_sizes" value="5,10,15,20,25,30,40,50">
        </div>
        <div class="field">
          <label for="cpm_thresholds">Peptide CPM cut-offs to test</label>
          <input type="text" id="cpm_thresholds" name="cpm_thresholds" value="10,15,20,30">
        </div>
        <div class="field">
          <label for="cpm_choice">CPM cut-off to use</label>
          <input type="text" id="cpm_choice" name="cpm_choice" placeholder="blank = first cut-off tested">
          <div class="hint">Used for mock-aware ranking and peptide presence.</div>
        </div>
        <div class="field">
          <label for="panel_choice">Number of top peptides to use</label>
          <input type="text" id="panel_choice" name="panel_choice" placeholder="blank = pipeline's tuned suggestion">
          <div class="hint">Leave blank to accept the tuned size the run works out for itself.</div>
        </div>
      </div>
      <div class="field">
        <label for="min_positive_peptides">Minimum positive peptides / hits</label>
        <input type="text" id="min_positive_peptides" name="min_positive_peptides" value="5">
        <div class="hint">Used for breadth-gated metrics, the default hybrid positive call, and
          heatmap sample suppression.</div>
      </div>
    </fieldset>

    <fieldset>
      <legend>Calling and display</legend>
      <div class="grid2">
        <div class="field">
          <label for="call_mode">Primary call</label>
          <select id="call_mode" name="call_mode"></select>
        </div>
        <div class="field">
          <label for="include_mocks">Mock columns in the heatmap</label>
          <select id="include_mocks" name="include_mocks">
            <option value="y">Include</option>
            <option value="n">Exclude</option>
          </select>
        </div>
        <div class="field" id="zclip_field" style="display:none">
          <label for="z_clip">Z-score display cutoff</label>
          <input type="text" id="z_clip" name="z_clip" value="3.5">
        </div>
        <div class="field">
          <label for="output_set">Output set</label>
          <select id="output_set" name="output_set">
            <option value="minimal">minimal — core outputs only</option>
            <option value="full">full — all ancillary outputs</option>
          </select>
        </div>
      </div>
      <div class="field">
        <label for="elisa_csv">ELISA CSV <span class="tag">optional</span></label>
        <div class="row">
          <input type="text" id="elisa_csv" name="elisa_csv" placeholder="sample name, 1/0">
          <button class="small" type="button" onclick="openBrowser('elisa_csv','csv')">Browse</button>
        </div>
      </div>
      <div class="field">
        <label for="extra_args">Additional pipeline arguments <span class="tag">optional</span></label>
        <input type="text" id="extra_args" name="extra_args" placeholder="--hybrid-min-positive-hits 3">
        <div class="hint">Must start with a flag. Quoted values with spaces are fine here.</div>
      </div>
    </fieldset>

    <fieldset>
      <legend>Run-versus-run comparison</legend>
      <div class="field">
        <label for="compare_enabled">Compare with a second run directory afterwards?</label>
        <select id="compare_enabled" name="compare_enabled" onchange="syncCompare()">
          <option value="n">No</option>
          <option value="y">Yes — run compare_runs_coverage.R</option>
        </select>
        <div class="hint">Same step the shell wrapper offers, and the same limits apply:
          compare_runs_coverage.R reads a hardcoded <code>africa2_and_s__selected_panel.txt</code>
          from each directory and needs one <code>heatmap_positive_values_*.csv</code> per method,
          so it only works on the ranking-sweep directories built that way. Its two figures are
          written to the <em>parent</em> of the output directory; the GUI picks them up from there.</div>
      </div>
      <div id="compare_fields" style="display:none">
        <div class="field">
          <label for="compare_dir">Second run directory</label>
          <div class="row">
            <input type="text" id="compare_dir" name="compare_dir">
            <button class="small" type="button" onclick="openBrowser('compare_dir','dir')">Browse</button>
          </div>
        </div>
        <div class="grid2">
          <div class="field">
            <label for="label1">Label for this run</label>
            <input type="text" id="label1" name="label1" placeholder="defaults to output directory name">
          </div>
          <div class="field">
            <label for="label2">Label for second run</label>
            <input type="text" id="label2" name="label2" placeholder="defaults to directory name">
          </div>
        </div>
      </div>
    </fieldset>

    <div id="formmsg"></div>
    <pre id="cmdbox" class="cmd" style="display:none"></pre>

    <div class="actions">
      <button class="primary" type="button" id="runbtn" onclick="startRun()">Run pipeline</button>
      <button type="button" onclick="preview()">Preview command</button>
      <button type="button" id="stopbtn" onclick="stopRun()" disabled>Stop</button>
      <label style="margin-left:auto;font-size:12px;color:var(--muted)">
        <input type="checkbox" id="overwrite" name="overwrite"> allow overwrite
      </label>
    </div>
  </form>
  </div>

  <!-- --------------------------------------------------- results column -->
  <div class="col">
    <div class="panel">
      <div class="tabs">
        <button class="on" onclick="showTab('log',this)">Log</button>
        <button onclick="showTab('plots',this)">Plots <span id="nplots" class="tag">0</span></button>
        <button onclick="showTab('files',this)">Files <span id="nfiles" class="tag">0</span></button>
      </div>
      <div class="tabbody">
        <div id="tab-log"><pre id="log">No run yet.</pre></div>
        <div id="tab-plots" style="display:none"><div id="gallery" class="gallery"></div></div>
        <div id="tab-files" style="display:none"><div id="filelist"></div><div id="tablepreview"></div></div>
      </div>
    </div>
  </div>
</main>

<div class="modal" id="browser">
  <div class="box">
    <header><h1 id="browsertitle">Choose</h1></header>
    <div class="list" id="browserlist"></div>
    <div class="foot">
      <button type="button" onclick="chooseCurrentDir()" id="choosedirbtn">Use this directory</button>
      <button type="button" style="margin-left:auto" onclick="closeBrowser()">Cancel</button>
    </div>
  </div>
</div>

<div id="lightbox" onclick="this.classList.remove('on')"><img id="lightimg" alt=""></div>

<script>
const HEATMAP_MODES = __HEATMAP_MODES__;
const COLLAPSE_MODES = __COLLAPSE_MODES__;
const CALL_MODES = __CALL_MODES__;
const DEFAULT_OUT_DIR = __DEFAULT_OUT_DIR__;

const FIELDS = ["matrix_csv","filter_term","output_prefix","out_dir","heatmap_mode","mock_regex",
  "true_negative_controls","analysis_mode","comparison_regex","pos_regex","neg_regex","pos_select",
  "neg_select","collapse","collapse_mode","collapse_floor","panel_sizes","cpm_thresholds",
  "cpm_choice","panel_choice","min_positive_peptides","call_mode","include_mocks","z_clip",
  "output_set","elisa_csv","extra_args","compare_enabled","compare_dir","label1","label2"];

let runId = null, logOffset = 0, timer = null, lastImages = "";

function el(id){ return document.getElementById(id); }

function fillSelect(id, entries){
  const sel = el(id);
  sel.innerHTML = "";
  for (const [value, label] of entries){
    const o = document.createElement("option");
    o.value = value; o.textContent = label;
    sel.appendChild(o);
  }
}

function getForm(){
  const data = {};
  for (const name of FIELDS){ const node = el(name); if (node) data[name] = node.value; }
  data.overwrite = el("overwrite").checked ? "y" : "n";
  return data;
}

function setForm(data){
  for (const name of FIELDS){
    const node = el(name);
    if (node && data[name] !== undefined && data[name] !== null) node.value = data[name];
  }
  syncMode(); syncCollapse(); syncCompare(); syncZclip();
}

function syncMode(){
  const supervised = el("analysis_mode").value === "supervised";
  el("sup_fields").style.display = supervised ? "" : "none";
  el("unsup_fields").style.display = supervised ? "none" : "";
}
function syncCollapse(){
  const on = el("collapse").value === "y";
  el("collapse_fields").style.display = on ? "" : "none";
  el("floor_field").style.display = (on && el("collapse_mode").value !== "raw_cpm") ? "" : "none";
}
function syncCompare(){
  el("compare_fields").style.display = el("compare_enabled").value === "y" ? "" : "none";
}
function syncZclip(){
  el("zclip_field").style.display = el("heatmap_mode").value === "zscore" ? "" : "none";
}

function message(kind, text){
  el("formmsg").innerHTML = text ? '<div class="msg ' + kind + '"></div>' : "";
  if (text) el("formmsg").firstChild.textContent = text;
}

async function api(path, options){
  const res = await fetch(path, options);
  const body = await res.json();
  if (!res.ok) throw new Error(body.error || ("HTTP " + res.status));
  return body;
}

async function preview(){
  message("", "");
  try {
    const body = await api("/api/preview", {method:"POST", body: JSON.stringify(getForm())});
    el("cmdbox").style.display = "";
    let text = body.command;
    if (body.notes && body.notes.length) text += "\n\n" + body.notes.map(n => "note: " + n).join("\n");
    if (body.existing && body.existing.length)
      text += "\n\nwarning: " + body.existing.length + " file(s) matching " +
              el("output_prefix").value + "_* already exist here and will be overwritten.";
    el("cmdbox").textContent = text;
  } catch (err){ el("cmdbox").style.display = "none"; message("err", err.message); }
}

async function startRun(){
  message("", "");
  el("runbtn").disabled = true;
  try {
    const body = await api("/api/run", {method:"POST", body: JSON.stringify(getForm())});
    runId = body.id; logOffset = 0; lastImages = "";
    el("log").textContent = "";
    el("gallery").innerHTML = ""; el("filelist").innerHTML = ""; el("tablepreview").innerHTML = "";
    el("stopbtn").disabled = false;
    el("cmdbox").style.display = "";
    el("cmdbox").textContent = body.command;
    showTab("log", document.querySelector(".tabs button"));
    poll();
  } catch (err){
    el("runbtn").disabled = false;
    if (err.message.indexOf("already exist") >= 0)
      message("warn", err.message + "\n\nTick “allow overwrite” to go ahead.");
    else message("err", err.message);
  }
}

async function stopRun(){
  if (!runId) return;
  await fetch("/api/stop?id=" + encodeURIComponent(runId), {method:"POST"});
}

function setPill(state, elapsed){
  const pill = el("pill");
  pill.className = "pill " + state;
  const secs = elapsed ? " · " + Math.round(elapsed) + "s" : "";
  pill.textContent = state + secs;
}

async function poll(){
  if (!runId) return;
  try {
    const body = await api("/api/status?id=" + encodeURIComponent(runId) + "&offset=" + logOffset);
    if (body.log){
      const box = el("log");
      const stuck = box.scrollTop + box.clientHeight >= box.scrollHeight - 40;
      box.textContent += body.log;
      if (stuck) box.scrollTop = box.scrollHeight;
    }
    logOffset = body.offset;
    setPill(body.state, body.elapsed);
    renderOutputs(body.outputs);
    if (body.state === "running" || body.state === "starting" || body.state === "stopping"){
      timer = setTimeout(poll, 900);
    } else {
      el("runbtn").disabled = false;
      el("stopbtn").disabled = true;
      if (body.error) message("err", body.error);
      else if (body.state === "finished"){
        const n = body.outputs.images.length;
        message("", "");
        if (n) showTab("plots", document.querySelectorAll(".tabs button")[1]);
      }
    }
  } catch (err){
    el("runbtn").disabled = false; el("stopbtn").disabled = true;
    message("err", err.message);
  }
}

function human(bytes){
  if (bytes < 1024) return bytes + " B";
  if (bytes < 1048576) return (bytes/1024).toFixed(0) + " KB";
  return (bytes/1048576).toFixed(1) + " MB";
}

function fileUrl(name, root){
  return "/file?id=" + encodeURIComponent(runId) + "&name=" + encodeURIComponent(name)
       + "&root=" + (root || 0);
}

function renderOutputs(outputs){
  const images = outputs.images, tables = outputs.tables.concat(outputs.others);
  el("nplots").textContent = images.length;
  el("nfiles").textContent = images.length + tables.length;

  const key = images.map(i => i.root + ":" + i.name + ":" + i.mtime).join("|");
  if (key !== lastImages){
    lastImages = key;
    const gallery = el("gallery");
    gallery.innerHTML = "";
    if (!images.length){
      gallery.innerHTML = '<div class="empty">No figures yet.</div>';
    } else {
      for (const img of images){
        const url = fileUrl(img.name, img.root) + "&v=" + Math.round(img.mtime);
        const card = document.createElement("div");
        card.className = "plot";
        const cap = document.createElement("div");
        cap.className = "cap";
        const nm = document.createElement("span"); nm.className = "name"; nm.textContent = img.name;
        const dl = document.createElement("a"); dl.href = url; dl.download = img.name; dl.textContent = "download";
        const sz = document.createElement("span"); sz.className = "sz"; sz.textContent = human(img.size);
        cap.appendChild(nm); cap.appendChild(sz); cap.appendChild(dl);
        const im = document.createElement("img");
        im.src = url; im.alt = img.name; im.loading = "lazy";
        im.onclick = () => { el("lightimg").src = url; el("lightbox").classList.add("on"); };
        card.appendChild(cap); card.appendChild(im);
        gallery.appendChild(card);
      }
    }
  }

  const list = el("filelist");
  const all = images.concat(tables);
  let html = '<table class="files">';
  for (const f of all){
    const clickable = /[.](csv|tsv|txt|log)$/i.test(f.name);
    html += "<tr><td>" + (clickable
        ? '<span class="fname" data-name="' + escapeAttr(f.name) + '" data-root="' + (f.root || 0) + '">'
          + escapeHtml(f.name) + "</span>"
        : escapeHtml(f.name))
      + '</td><td class="sz">' + human(f.size) + "</td>"
      + '<td><a href="' + fileUrl(f.name, f.root) + '" download>download</a></td></tr>';
  }
  html += "</table>";
  list.innerHTML = all.length ? html : '<div class="empty">No output files yet.</div>';
  for (const node of list.querySelectorAll(".fname"))
    node.onclick = () => showTable(node.dataset.name, node.dataset.root);
}

function escapeHtml(s){ const d = document.createElement("div"); d.textContent = s; return d.innerHTML; }
function escapeAttr(s){ return escapeHtml(s).replace(/"/g, "&quot;"); }

async function showTable(name, root){
  try {
    const body = await api("/api/table?id=" + encodeURIComponent(runId)
      + "&name=" + encodeURIComponent(name) + "&root=" + (root || 0));
    let html = "<h3 style='margin:16px 0 0;font-size:13px'>" + escapeHtml(name) + "</h3>";
    html += '<div class="tablewrap"><table>';
    if (body.header){
      html += "<tr>" + body.header.map(h => "<th>" + escapeHtml(h) + "</th>").join("") + "</tr>";
    }
    for (const row of body.rows)
      html += "<tr>" + row.map(c => "<td>" + escapeHtml(c) + "</td>").join("") + "</tr>";
    html += "</table></div>";
    if (body.truncated) html += '<div class="hint">First rows only — download for the full file.</div>';
    el("tablepreview").innerHTML = html;
  } catch (err){ el("tablepreview").innerHTML = '<div class="msg err">' + escapeHtml(err.message) + "</div>"; }
}

function showTab(name, button){
  for (const id of ["log","plots","files"]) el("tab-" + id).style.display = (id === name) ? "" : "none";
  for (const b of document.querySelectorAll(".tabs button")) b.classList.remove("on");
  if (button) button.classList.add("on");
}

// ------------------------------------------------------- sample columns
async function loadColumns(){
  el("colnote").textContent = "reading header...";
  el("samplebox").innerHTML = "";
  try {
    const params = new URLSearchParams({path: el("matrix_csv").value, mock_regex: el("mock_regex").value});
    const body = await api("/api/columns?" + params.toString());
    el("colnote").textContent = body.n_samples + " sample columns · " +
      body.mocks.length + " mock · " + body.nonmocks.length + " non-mock";
    let html = '<div class="samplebox"><table>';
    for (const s of body.nonmocks){
      html += '<tr><td class="i">' + s.index + "</td><td>" + escapeHtml(s.name) + "</td>"
        + '<td class="b">'
        + '<button type="button" class="small" data-i="' + s.index + '" data-t="pos_select">+pos</button> '
        + '<button type="button" class="small" data-i="' + s.index + '" data-t="neg_select">+neg</button>'
        + "</td></tr>";
    }
    for (const m of body.mocks)
      html += '<tr class="mockrow"><td class="i">—</td><td>' + escapeHtml(m) + "</td><td class='b'>mock</td></tr>";
    html += "</table></div>";
    el("samplebox").innerHTML = html;
    for (const b of el("samplebox").querySelectorAll("button"))
      b.onclick = () => addIndex(b.dataset.t, b.dataset.i);
  } catch (err){ el("colnote").textContent = ""; el("samplebox").innerHTML =
      '<div class="msg err" style="margin:8px 0">' + escapeHtml(err.message) + "</div>"; }
}

function addIndex(target, index){
  const node = el(target);
  const parts = node.value.split(",").map(s => s.trim()).filter(Boolean);
  if (!parts.includes(index)) parts.push(index);
  node.value = parts.join(",");
  if (el("analysis_mode").value !== "supervised"){
    el("analysis_mode").value = "supervised"; syncMode();
  }
}

// ------------------------------------------------------- file browser
let browserTarget = null, browserWant = "csv", browserDir = DEFAULT_OUT_DIR;

async function openBrowser(target, want){
  browserTarget = target; browserWant = want;
  const current = el(target).value;
  browserDir = current || DEFAULT_OUT_DIR;
  el("browsertitle").textContent = want === "dir" ? "Choose a directory" : "Choose a file";
  el("choosedirbtn").style.display = want === "dir" ? "" : "none";
  el("browser").classList.add("on");
  await listDir(browserDir, true);
}

function closeBrowser(){ el("browser").classList.remove("on"); }

async function listDir(path, guess){
  try {
    const params = new URLSearchParams({dir: path, want: browserWant, guess: guess ? "1" : "0"});
    const body = await api("/api/browse?" + params.toString());
    browserDir = body.dir;
    let html = '<div data-d="' + escapeAttr(body.parent) + '"><strong>..</strong> ' + escapeHtml(body.dir) + "</div>";
    for (const d of body.dirs)
      html += '<div data-d="' + escapeAttr(body.dir.replace(/\/$/, "") + "/" + d) + '">\u{1F4C1} ' + escapeHtml(d) + "</div>";
    for (const f of body.files)
      html += '<div data-f="' + escapeAttr(body.dir.replace(/\/$/, "") + "/" + f) + '">\u{1F4C4} ' + escapeHtml(f) + "</div>";
    el("browserlist").innerHTML = html;
    for (const node of el("browserlist").querySelectorAll("div")){
      node.onclick = () => {
        if (node.dataset.d) listDir(node.dataset.d, false);
        else { el(browserTarget).value = node.dataset.f; closeBrowser(); }
      };
    }
  } catch (err){ el("browserlist").innerHTML = '<div class="msg err">' + escapeHtml(err.message) + "</div>"; }
}

function chooseCurrentDir(){ el(browserTarget).value = browserDir; closeBrowser(); }

// ------------------------------------------------------------- start-up
(function init(){
  fillSelect("heatmap_mode", HEATMAP_MODES);
  fillSelect("collapse_mode", COLLAPSE_MODES.map(m => [m, m]));
  fillSelect("call_mode", CALL_MODES.map(m => [m, m]));
  el("heatmap_mode").value = "mock_q75";
  el("out_dir").value = DEFAULT_OUT_DIR;
  el("heatmap_mode").addEventListener("change", syncZclip);
  fetch("/api/last").then(r => r.json()).then(data => { if (data && data.matrix_csv) setForm(data); })
    .catch(() => {});
  syncMode(); syncCollapse(); syncCompare(); syncZclip();
})();
</script>
</body>
</html>
"""


def format_command(cmd):
    """One line per flag, keeping each flag with its value, for the preview."""
    # cmd always opens with interpreter, -u, the pipeline and the matrix; those
    # belong together on the first line.
    lines = [" ".join(shlex.quote(c) for c in cmd[:4])]
    i = 4
    while i < len(cmd):
        token = cmd[i]
        if token.startswith("-") and i + 1 < len(cmd) and not cmd[i + 1].startswith("-"):
            lines.append("%s %s" % (shlex.quote(token), shlex.quote(cmd[i + 1])))
            i += 2
        else:
            lines.append(shlex.quote(token))
            i += 1
    return " \\\n    ".join(lines)


def render_page():
    return (PAGE
            .replace("__HEATMAP_MODES__", json.dumps([[k, v] for k, v in HEATMAP_MODES.items()]))
            .replace("__COLLAPSE_MODES__", json.dumps(COLLAPSE_MODES))
            .replace("__CALL_MODES__", json.dumps(CALL_MODES))
            .replace("__DEFAULT_OUT_DIR__", json.dumps(HERE)))


# ===========================================================================
# HTTP
# ===========================================================================
class Handler(BaseHTTPRequestHandler):
    server_version = "PhIPseqGUI/1.0"

    def log_message(self, fmt, *args):
        if VERBOSE:
            sys.stderr.write("%s - %s\n" % (self.address_string(), fmt % args))

    # -- plumbing ----------------------------------------------------------
    def _send(self, code, body, ctype="application/json; charset=utf-8", extra=None):
        if isinstance(body, str):
            body = body.encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        for key, value in (extra or {}).items():
            self.send_header(key, value)
        self.end_headers()
        self.wfile.write(body)

    def _json(self, payload, code=200):
        self._send(code, json.dumps(payload))

    def _error(self, message, code=400):
        self._json({"error": str(message)}, code=code)

    def _body(self):
        length = int(self.headers.get("Content-Length") or 0)
        raw = self.rfile.read(length) if length else b"{}"
        try:
            return json.loads(raw.decode("utf-8") or "{}")
        except ValueError:
            raise FormError("Could not read the submitted form.")

    def _run_from_query(self, query):
        run = RUNS.get((query.get("id") or [""])[0])
        if run is None:
            raise FormError("Unknown run id — the GUI may have been restarted.")
        return run

    @staticmethod
    def _resolve_output(run, name, root=0):
        """Resolve a requested filename strictly inside one of the run's roots."""
        try:
            base = os.path.realpath(run.roots[int(root)])
        except (ValueError, IndexError):
            raise FormError("Unknown output location.")
        target = os.path.realpath(os.path.join(base, name))
        if target != base and not target.startswith(base + os.sep):
            raise FormError("Refusing to serve a file outside the output directory.")
        if not os.path.isfile(target):
            raise FormError("File not found: %s" % name)
        return target

    # -- routes ------------------------------------------------------------
    def do_GET(self):
        parsed = urlparse(self.path)
        route, query = parsed.path, parse_qs(parsed.query)
        try:
            if route == "/":
                self._send(200, render_page(), "text/html; charset=utf-8")
            elif route == "/api/last":
                self._json(load_last_form())
            elif route == "/api/columns":
                self._json(read_sample_columns((query.get("path") or [""])[0],
                                               (query.get("mock_regex") or ["mock"])[0]))
            elif route == "/api/browse":
                path = (query.get("dir") or [HERE])[0]
                want = (query.get("want") or ["csv"])[0]
                # A path typed into a file field is a file, not a directory;
                # open its containing directory rather than falling back home.
                if (query.get("guess") or ["0"])[0] == "1" and path and not os.path.isdir(_abspath(path)):
                    path = os.path.dirname(_abspath(path))
                self._json(browse_dir(path, want))
            elif route == "/api/status":
                run = self._run_from_query(query)
                offset = int((query.get("offset") or ["0"])[0])
                self._json(run.status(offset))
            elif route == "/api/table":
                run = self._run_from_query(query)
                path = self._resolve_output(run, (query.get("name") or [""])[0],
                                            (query.get("root") or ["0"])[0])
                self._json(preview_table(path))
            elif route == "/file":
                run = self._run_from_query(query)
                path = self._resolve_output(run, (query.get("name") or [""])[0],
                                            (query.get("root") or ["0"])[0])
                ctype = mimetypes.guess_type(path)[0] or "application/octet-stream"
                with open(path, "rb") as fh:
                    self._send(200, fh.read(), ctype)
            else:
                self._error("Not found", 404)
        except FormError as exc:
            self._error(exc)
        except Exception as exc:  # keep the server up whatever the browser asks for
            self._error("%s: %s" % (type(exc).__name__, exc), 500)

    def do_POST(self):
        parsed = urlparse(self.path)
        route, query = parsed.path, parse_qs(parsed.query)
        try:
            if route == "/api/preview":
                form = self._body()
                cmd, stdin_answers, out_dir, notes = build_command(form)
                self._json({
                    "command": "cd %s\n%s" % (shlex.quote(out_dir), format_command(cmd)),
                    "notes": notes,
                    "existing": existing_prefix_files(out_dir, _s(form, "output_prefix")),
                })
            elif route == "/api/run":
                form = self._body()
                cmd, stdin_answers, out_dir, notes = build_command(form)
                existing = existing_prefix_files(out_dir, _s(form, "output_prefix"))
                if existing and _s(form, "overwrite") not in ("y", "yes", "on", "true", "1"):
                    raise FormError(
                        "%d file(s) matching %s_* already exist in %s and would be overwritten."
                        % (len(existing), _s(form, "output_prefix"), out_dir))
                save_last_form(form)
                run_id = uuid.uuid4().hex[:12]
                run = Run(run_id, cmd, stdin_answers, out_dir, form, notes)
                with RUNS_LOCK:
                    RUNS[run_id] = run
                run.start()
                self._json({"id": run_id, "out_dir": out_dir,
                            "command": " ".join(shlex.quote(c) for c in cmd),
                            "notes": notes})
            elif route == "/api/stop":
                run = self._run_from_query(query)
                self._json({"stopped": run.stop()})
            else:
                self._error("Not found", 404)
        except FormError as exc:
            self._error(exc)
        except Exception as exc:
            self._error("%s: %s" % (type(exc).__name__, exc), 500)


VERBOSE = False


def main():
    global PYTHON_BIN, VERBOSE
    parser = argparse.ArgumentParser(
        description="Browser GUI for the PhIP-seq rank/call/plot pipeline v1.7.2.")
    parser.add_argument("--host", default="127.0.0.1",
                        help="Interface to bind. Default 127.0.0.1 (this machine only).")
    parser.add_argument("--port", type=int, default=8765, help="Port to listen on. Default 8765.")
    parser.add_argument("--python", default="python3",
                        help="Interpreter used to run the pipeline. Default python3.")
    parser.add_argument("--open", action="store_true", help="Open a browser window on start.")
    parser.add_argument("--verbose", action="store_true", help="Log every HTTP request.")
    args = parser.parse_args()

    PYTHON_BIN = args.python
    VERBOSE = args.verbose

    if not os.path.isfile(PY_SCRIPT):
        sys.exit("ERROR: Cannot find the pipeline at: %s" % PY_SCRIPT)
    if not shutil.which(PYTHON_BIN):
        sys.exit("ERROR: Interpreter not found on PATH: %s" % PYTHON_BIN)

    server = ThreadingHTTPServer((args.host, args.port), Handler)
    url = "http://%s:%d/" % ("127.0.0.1" if args.host in ("0.0.0.0", "") else args.host, args.port)
    print("PhIP-seq GUI ready")
    print("  pipeline : %s" % PY_SCRIPT)
    print("  python   : %s" % shutil.which(PYTHON_BIN))
    print("  open     : %s" % url)
    if args.host not in ("127.0.0.1", "localhost"):
        print("  NOTE: bound to %s — anything that can reach this port can start runs "
              "and read files under the output directory." % args.host)
    print("  (over SSH: ssh -L %d:127.0.0.1:%d <this-host>, then open the URL locally)"
          % (args.port, args.port))
    print("Ctrl-C to stop.")
    sys.stdout.flush()   # so the URL still shows when the banner is piped to a log
    if args.open:
        import webbrowser
        threading.Thread(target=lambda: (time.sleep(0.6), webbrowser.open(url)), daemon=True).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down.")
        for run in list(RUNS.values()):
            run.stop()
        server.shutdown()


if __name__ == "__main__":
    main()
