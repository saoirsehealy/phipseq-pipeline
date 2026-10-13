#!/usr/bin/env python3
# =============================================================================
# PhIP-seq virus preparation + rank/call GUI  (remote / SSH edition)
# -----------------------------------------------------------------------------
# A browser front-end that logs in to a remote Linux server with the user's own
# account, then chains the two halves of the workflow on that server:
#
#   1. annotate_plates.py          phipflow wide_data -> real sample and peptide
#                                  names, plates combined   (optional)
#   2. prepare_virus_matrix_v2.py  library -> virus peptides -> BLAST ->
#                                  coordinates -> genotypes -> relabelled matrix
#   3. phipseq_rank_call_and_plot  rank -> tune panel -> call -> heatmap
#   4. plot_genome.py              where the signal sits along the genome
#
# Each stage feeds the next automatically: the combined matrix into stage two,
# the relabelled matrix into stage three, and that matrix plus the panel stage
# three selected into stage four. Every command is shown in full before anything
# runs. Log output streams back live and the figures and tables produced are
# listed and viewable in the browser.
#
# Two constraints the form enforces, because they fail confusingly otherwise:
#
#   * Stage two takes EITHER --search-term (a library that names your virus
#     reliably, such as OCCUSCAN) OR --sequence-search (one that does not, such
#     as RESPSCAN). v2 refuses both together.
#   * Stage four needs a panel from an UNCOLLAPSED stage three. Collapsing
#     replaces peptide names with bare positions, which carry no protein or
#     coordinates, so the genome plot has nothing to place.
#
# CREDENTIALS
# -----------
# The username and password are the person's ordinary Linux login for the
# server they name. They are used once, to open an SSH connection, and are then
# discarded: only the live SSH client object is kept, keyed to an opaque session
# cookie. The password is never written to disk, never logged, and never sent
# back to the browser.
#
#   *** The browser -> GUI hop is plain HTTP. That hop would carry a real
#   *** account password, so by default this server BINDS TO LOOPBACK and
#   *** REFUSES credentials from any other address. Reach it through an SSH
#   *** tunnel:
#   ***     ssh -L 8770:127.0.0.1:8770 <host-running-this>
#   *** --allow-insecure-credentials lifts the refusal. Do not use it on an
#   *** untrusted network, and never on the public internet: anyone in path
#   *** could read the password.
#
# Host keys are verified against ~/.ssh/known_hosts and the system host keys,
# and an unknown host is refused rather than silently trusted. Connect once
# with `ssh user@host` to record the key, or pass --auto-add-host-key if you
# accept the first-use risk.
#
# Every remote path is normalised on the server (so symlinks and ".." cannot
# escape) and rejected unless it sits inside the root directory chosen at login.
#
# Usage:
#   python3 phipseq_prepare_gui.py
#   python3 phipseq_prepare_gui.py --port 8770 --open
#   python3 phipseq_prepare_gui.py --host 0.0.0.0 --allow-insecure-credentials
# =============================================================================
import argparse
import html
import json
import mimetypes
import os
import posixpath
import re
import secrets
import shlex
import socket
import sys
import threading
import time
import webbrowser
from http import HTTPStatus
from http.cookies import SimpleCookie
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

try:
    import paramiko
except ImportError:
    raise SystemExit("ERROR: paramiko is required.\n"
                     "       python3 -m pip install --user paramiko")

SESSION_IDLE_SECONDS = 4 * 3600
VIEWABLE = {".png", ".csv", ".tsv", ".txt", ".log", ".fa", ".fasta", ".pdf"}
HEATMAP_MODES = ["raw_cpm", "zscore", "minus_mean", "positive_mock_median",
                 "mock_q75", "mock_q100", "mock_q100_plus_median", "mean_3_5sd"]
# The login form is prefilled from these. Nothing here is binding -- every one
# is an editable field -- but a site that sets the environment variables below
# (in /etc/profile.d, a module file, or the launcher) gives its own users
# sensible defaults instead of whoever set this copy up first.
#
#   PHIPSEQ_PIPELINE_DIR   folder holding prepare_virus_matrix_v2.py et al.
#   PHIPSEQ_ROOT           project area the session is confined to
#   PHIPSEQ_PYTHON         interpreter to run the scripts with on the server
#
# PHIPSEQ_PIPELINE_DIR falls back to this file's own directory, which is right
# whenever the GUI is kept beside the scripts it drives -- the normal case.
DEFAULT_PIPELINE_DIR = os.environ.get(
    "PHIPSEQ_PIPELINE_DIR", str(Path(__file__).resolve().parent))
DEFAULT_ROOT = os.environ.get("PHIPSEQ_ROOT", DEFAULT_PIPELINE_DIR)
DEFAULT_PYTHON = os.environ.get("PHIPSEQ_PYTHON", "python3")
# Blank by default. A site's own analysis server goes in PHIPSEQ_SERVER (set in
# the launcher), rather than one institution's hostname being the default every
# other institution has to delete.
DEFAULT_SERVER = os.environ.get("PHIPSEQ_SERVER", "")
# Public git URL the scripts are installed from. With this set, a colleague
# needs nothing but this one file: the GUI clones the pipeline onto their own
# server over the SSH connection they just opened, so they never touch git.
DEFAULT_REPO = os.environ.get(
    "PHIPSEQ_REPO", "https://github.com/saoirsehealy/phipseq-pipeline.git")

ALLOW_INSECURE = False
AUTO_ADD_HOST_KEY = False


class FormError(Exception):
    pass


# ------------------------------------------------------------------- session
class Session:
    """One authenticated SSH connection plus the state of its current run."""

    def __init__(self, client, host, user, root, pipeline_dir, python_exe):
        self.client = client
        self.host = host
        self.user = user
        self.root = root
        self.pipeline_dir = pipeline_dir
        self.python = python_exe
        self.repo = ""
        self.seen = time.monotonic()
        self._sftp = None
        self._sftp_lock = threading.Lock()
        # run state
        self.lock = threading.Lock()
        self.lines = []
        self.running = False
        self.done = False
        self.failed = False
        self.out_dir = None

    def sftp(self):
        with self._sftp_lock:
            if self._sftp is None:
                self._sftp = self.client.open_sftp()
            return self._sftp

    def close(self):
        try:
            if self._sftp is not None:
                self._sftp.close()
        except Exception:
            pass
        try:
            self.client.close()
        except Exception:
            pass

    # -- run state
    def reset(self, out_dir):
        with self.lock:
            self.lines, self.running, self.done, self.failed = [], True, False, False
            self.out_dir = out_dir

    def add(self, text):
        with self.lock:
            for line in str(text).splitlines():
                self.lines.append(line)

    def finish(self, failed):
        with self.lock:
            self.running, self.done, self.failed = False, True, failed

    def read(self, since):
        with self.lock:
            return self.lines[since:], len(self.lines), self.running, self.done, self.failed


SESSIONS = {}
SESSIONS_LOCK = threading.Lock()
FAILURES = {}
FAILURES_LOCK = threading.Lock()


def get_session(token):
    if not token:
        return None
    now = time.monotonic()
    with SESSIONS_LOCK:
        sess = SESSIONS.get(token)
        if sess is None:
            return None
        if now - sess.seen > SESSION_IDLE_SECONDS:
            SESSIONS.pop(token, None)
            sess.close()
            return None
        sess.seen = now
        return sess


def drop_session(token):
    with SESSIONS_LOCK:
        sess = SESSIONS.pop(token, None)
    if sess:
        sess.close()


def throttle(addr):
    with FAILURES_LOCK:
        count, until = FAILURES.get(addr, (0, 0.0))
    return max(0.0, until - time.monotonic())


def note_failure(addr):
    with FAILURES_LOCK:
        count, _ = FAILURES.get(addr, (0, 0.0))
        count += 1
        delay = 0.0 if count < 3 else min(300.0, 2.0 ** (count - 2))
        FAILURES[addr] = (count, time.monotonic() + delay)


def clear_failures(addr):
    with FAILURES_LOCK:
        FAILURES.pop(addr, None)


# ---------------------------------------------------------------- ssh helpers
def ssh_connect(host, port, user, password):
    client = paramiko.SSHClient()
    client.load_system_host_keys()
    known = Path.home() / ".ssh" / "known_hosts"
    if known.exists():
        try:
            client.load_host_keys(str(known))
        except Exception:
            pass
    client.set_missing_host_key_policy(
        paramiko.AutoAddPolicy() if AUTO_ADD_HOST_KEY else paramiko.RejectPolicy())
    client.connect(hostname=host, port=port, username=user, password=password,
                   look_for_keys=False, allow_agent=False, timeout=20,
                   auth_timeout=20, banner_timeout=20)
    return client


REQUIRED_SCRIPTS = ("prepare_virus_matrix_v2.py",
                    "phipseq_rank_call_and_plot_v1.7.2.py",
                    "plot_genome.py")


def check_install(sess):
    """What is missing on the server, as (scripts, packages, binaries) lists.

    Called at login and again after every repair, so the setup page always
    reflects the server as it is now rather than as it was when we logged in.
    """
    scripts = [name for name in REQUIRED_SCRIPTS
               if remote_run(sess, "test -f " + shlex.quote(
                   posixpath.join(sess.pipeline_dir, name)))[0] != 0]

    packages = []
    st, out = remote_run(sess, shlex.quote(sess.python) + " -c " + shlex.quote(
        "import pandas, numpy, matplotlib, scipy"))
    if st != 0:
        found = re.findall(r"No module named '([^']+)'", out or "")
        packages = found or [f"cannot run {sess.python}: "
                             f"{(out or 'not found').strip().splitlines()[-1][:160]}"]

    binaries = [tool for tool in ("blastp", "tblastn", "makeblastdb")
                if remote_run(sess, "command -v " + tool)[0] != 0]
    return scripts, packages, binaries


def remote_run(sess, command, timeout=60):
    """Run a short command, return (exit_status, combined output)."""
    chan = sess.client.get_transport().open_session(timeout=timeout)
    chan.settimeout(timeout)
    chan.set_combine_stderr(True)
    chan.exec_command(command)
    out = []
    try:
        while True:
            data = chan.recv(65536)
            if not data:
                break
            out.append(data.decode("utf-8", "replace"))
    except socket.timeout:
        out.append("\n[timed out]")
    status = chan.recv_exit_status()
    chan.close()
    return status, "".join(out)


def safe_remote(sess, raw, must_exist=True, want_dir=False):
    """Normalise a remote path ON THE SERVER and refuse anything outside root."""
    text = (raw or "").strip()
    if not text:
        raise FormError("A required path was left blank.")
    candidate = text if text.startswith("/") else posixpath.join(sess.root, text)
    candidate = posixpath.normpath(candidate)
    sftp = sess.sftp()
    # normalize() resolves on the server, so symlinks cannot be used to escape.
    # An output directory may not exist yet, and neither may several levels of its
    # parents, so walk up to the deepest ancestor that does exist, resolve that,
    # then re-append the components below it.
    trailing = []
    probe = candidate
    while True:
        try:
            resolved = posixpath.join(sftp.normalize(probe), *reversed(trailing)) \
                if trailing else sftp.normalize(probe)
            break
        except IOError:
            parent = posixpath.dirname(probe)
            if not parent or parent == probe:
                raise FormError(f"Cannot resolve path: {text}")
            trailing.append(posixpath.basename(probe))
            probe = parent
            if len(trailing) > 32:
                raise FormError(f"Path is too deeply nested to resolve: {text}")
    root = sess.root.rstrip("/")
    if resolved != root and not resolved.startswith(root + "/"):
        raise FormError(f"Path is outside the permitted root ({root}): {text}")
    if must_exist:
        try:
            attrs = sftp.stat(resolved)
        except IOError:
            raise FormError(f"Does not exist on {sess.host}: {text}")
        import stat as statmod
        is_dir = statmod.S_ISDIR(attrs.st_mode)
        if want_dir and not is_dir:
            raise FormError(f"Not a directory: {text}")
        if not want_dir and is_dir:
            raise FormError(f"Not a file: {text}")
    return resolved


def remote_listing(sess, path):
    import stat as statmod
    sftp = sess.sftp()
    entries = []
    for attr in sftp.listdir_attr(path):
        if attr.filename.startswith("."):
            continue
        entries.append({"name": attr.filename,
                        "dir": statmod.S_ISDIR(attr.st_mode),
                        "size": attr.st_size})
    entries.sort(key=lambda e: (not e["dir"], e["name"].lower()))
    return entries


# ------------------------------------------------------------ command building
def build_commands(sess, form):
    virus = (form.get("virus") or "").strip()
    if not re.fullmatch(r"[A-Za-z0-9._-]+", virus or ""):
        raise FormError("Virus label is required (letters, digits, . _ - only).")
    # v2 offers two ways to find the virus's peptides and refuses both at once.
    search = (form.get("search_term") or "").strip()
    sequence_search = bool(form.get("sequence_search"))
    if sequence_search and search:
        raise FormError("Choose EITHER a search term OR sequence search, not both.")
    if not sequence_search and not search:
        raise FormError("Give a search term, or tick sequence search for a library "
                        "whose peptide names cannot be trusted.")

    prepare = posixpath.join(sess.pipeline_dir, "prepare_virus_matrix_v2.py")
    annotate = posixpath.join(sess.pipeline_dir, "annotate_plates.py")
    pipeline = posixpath.join(sess.pipeline_dir,
                              "phipseq_rank_call_and_plot_v1.7.2.py")

    coord_ref = safe_remote(sess, form.get("coord_ref"))
    out_dir = safe_remote(sess, form.get("out_dir"), must_exist=False)

    steps = []
    if form.get("run_annotate"):
        raw_plates = [line.strip() for line in
                      (form.get("plate_dirs") or "").replace(",", "\n").splitlines()
                      if line.strip()]
        if not raw_plates:
            raise FormError("Give at least one phipflow results folder, one per line.")
        plate_paths = [safe_remote(sess, plate, want_dir=True) for plate in raw_plates]
        if len(set(plate_paths)) != len(plate_paths):
            raise FormError("The same plate folder was listed more than once.")
        combined_dir = safe_remote(sess, form.get("combined_dir") or
                                   posixpath.join(out_dir, "combined"), must_exist=False)
        steps.append(("Annotate and combine plates",
                      [sess.python, "-u", annotate] + plate_paths
                      + ["--out-dir", combined_dir],
                      combined_dir))
        # Stage one supplies the combined CPM matrix with `peptide` as its first
        # column, which is what v2 expects. v2 reads peptide SEQUENCES only from
        # the translated library, so that file is required on both paths.
        cpm = posixpath.join(combined_dir, "data_cpm_final_annotated.csv")
    else:
        cpm = safe_remote(sess, form.get("cpm"))
    sequence_source = ["--library", safe_remote(sess, form.get("library"))]

    segment = (form.get("segment") or "").strip()
    if segment and not re.fullmatch(r"[A-Za-z0-9._-]+", segment):
        raise FormError("Segment may use letters, digits, . _ - only.")
    prefix = virus + (f"_{segment}" if segment else "")
    # v2 names the output file, not a directory. Keep a .csv extension: the
    # extension decides the delimiter, and the rank/call step wants a CSV.
    matrix = posixpath.join(out_dir, f"{prefix}_relabelled.csv")

    prep = ([sess.python, "-u", prepare, "--virus", virus] + sequence_source
            + ["--coord-ref", coord_ref, "--cpm", cpm, "--output", matrix])
    prep += ["--sequence-search"] if sequence_search else ["--search-term", search]
    if segment:
        prep += ["--segment", segment]
    for field_name, flag in (("coord_map", "--coord-map"), ("geno_ref", "--geno-ref")):
        if (form.get(field_name) or "").strip():
            prep += [flag, safe_remote(sess, form.get(field_name))]
    if (form.get("cpm_name_column") or "").strip():
        prep += ["--cpm-name-column", form["cpm_name_column"].strip()]
    if (form.get("cpm_drop_columns") or "").strip():
        prep += ["--cpm-drop-columns", form["cpm_drop_columns"].strip()]
    # v2 keeps unselected rows by default (so the matrix stays usable for other
    # viruses); the old --keep-all-rows is now the absence of --drop-unselected.
    if form.get("drop_unselected"):
        prep += ["--drop-unselected"]
    for field_name, flag, kind in (("min_query_coverage", "--min-query-coverage", float),
                                   ("evalue", "--evalue", str),
                                   ("min_identity", "--min-identity", float),
                                   ("ambiguity_fraction", "--ambiguity-fraction", float)):
        raw = (form.get(field_name) or "").strip()
        if raw:
            try:
                kind(raw)
            except ValueError:
                raise FormError(f"{flag} must be a {kind.__name__}.")
            prep += [flag, raw]

    steps.append(("Prepare virus matrix", prep, out_dir))
    if not form.get("run_pipeline"):
        return steps, out_dir, matrix

    run_dir = safe_remote(sess, form.get("run_dir") or posixpath.join(out_dir, "run"),
                          must_exist=False)
    mode = (form.get("heatmap_mode") or "mock_q75").strip()
    if mode not in HEATMAP_MODES:
        raise FormError(f"Unknown heatmap/ranking mode: {mode}")
    pipe = [sess.python, "-u", pipeline, matrix,
            "--filter-term", (form.get("filter_term") or virus).strip(),
            "--output-prefix", (form.get("output_prefix") or prefix).strip(),
            "--mock-regex", (form.get("mock_regex") or "mock").strip(),
            "--heatmap-mode", mode]
    pipe += ["--supervised"] if form.get("supervised") else ["--unsupervised"]
    for field_name, flag in (("pos_regex", "--pos-regex"), ("neg_regex", "--neg-regex"),
                             ("panel_sizes", "--panel-sizes"),
                             ("cpm_thresholds", "--cpm-thresholds")):
        if (form.get(field_name) or "").strip():
            pipe += [flag, form[field_name].strip()]
    for field_name, flag, kind in (("panel_choice", "--panel-choice", int),
                                   ("cpm_choice", "--cpm-choice", float)):
        raw = (form.get(field_name) or "").strip()
        if raw:
            try:
                kind(raw)
            except ValueError:
                raise FormError(f"{flag} must be a {kind.__name__}.")
            pipe += [flag, raw]
    pipe += (["--collapse-by-location-max"] if form.get("collapse")
             else ["--no-collapse-by-location-max"])
    extra = (form.get("extra_args") or "").strip()
    if extra:
        tokens = shlex.split(extra)
        if tokens and not tokens[0].startswith("-"):
            raise FormError(f"'{tokens[0]}' is not a flag — pipeline arguments start with '-'.")
        pipe += tokens

    steps.append(("Rank, call and plot", pipe, run_dir))

    if not form.get("run_genome"):
        return steps, run_dir, matrix

    # plot_genome.py reads the panel the rank/call step just wrote. A panel from
    # a COLLAPSED run holds bare position labels with no ORF or coordinates, so
    # the plot cannot use it.
    if form.get("collapse"):
        raise FormError("Genome plots need a panel from an uncollapsed run. "
                        "Untick location collapse, or untick genome plots.")
    out_prefix = (form.get("output_prefix") or prefix).strip()
    genome = posixpath.join(sess.pipeline_dir, "plot_genome.py")
    plots_dir = safe_remote(sess, form.get("genome_dir")
                            or posixpath.join(run_dir, "genome_plots"),
                            must_exist=False)
    measures = [m for m in (form.get("measures") or "pct zscore").split() if m]
    bad = [m for m in measures if m not in ("pct", "zscore", "cpm", "q75")]
    if bad:
        raise FormError(f"Unknown --measure value(s): {', '.join(bad)}")

    groups = [line.strip() for line in
              (form.get("groups") or "").replace(",", "\n").splitlines() if line.strip()]
    if not groups:
        raise FormError("Genome plots need at least one group, as Label=REGEX.")
    for spec in groups:
        if "=" not in spec or not spec.split("=", 1)[0].strip():
            raise FormError(f"Group {spec!r} must look like Label=REGEX.")

    plot = [sess.python, "-u", genome,
            "--matrix", matrix,
            "--panel", posixpath.join(run_dir, f"{out_prefix}_selected_panel.txt"),
            "--out-dir", plots_dir,
            "--virus", virus + (f" {segment}" if segment else ""),
            "--coordinate-level", "nt" if form.get("coord_map") else "aa",
            "--mock-regex", (form.get("mock_regex") or "mock").strip(),
            "--measure"] + measures
    for spec in groups:
        plot += ["--group", spec]
    for field_name, flag, kind in (("genome_length", "--genome-length", float),
                                   ("z_threshold", "--z-threshold", float)):
        raw = (form.get(field_name) or "").strip()
        if raw:
            try:
                kind(raw)
            except ValueError:
                raise FormError(f"{flag} must be a {kind.__name__}.")
            plot += [flag, raw]
    if (form.get("genome_label") or "").strip():
        plot += ["--genome-label", form["genome_label"].strip()]

    steps.append(("Genome plots", plot, plots_dir))
    return steps, plots_dir, matrix


# ------------------------------------------------------------------- the run
def execute(sess, steps, matrix):
    sess.reset(steps[-1][2])
    failed = False
    try:
        for index, (title, cmd, workdir) in enumerate(steps, start=1):
            sess.add(f"\n=== Step {index}/{len(steps)}: {title}   [on {sess.host}]")
            sess.add("  " + " ".join(shlex.quote(part) for part in cmd))
            sess.add("")
            status, _ = remote_run(sess, "mkdir -p " + shlex.quote(workdir))
            if status != 0:
                sess.add(f"ERROR: could not create {workdir}")
                failed = True
                break
            if title.startswith("Rank"):
                st, _ = remote_run(sess, "test -f " + shlex.quote(matrix))
                if st != 0:
                    sess.add(f"ERROR: expected matrix not found: {matrix}")
                    failed = True
                    break
            shell = ("cd " + shlex.quote(workdir) + " && "
                     + " ".join(shlex.quote(part) for part in cmd))
            chan = sess.client.get_transport().open_session()
            chan.set_combine_stderr(True)
            chan.exec_command(shell)
            # The pipeline asks a few questions with no command-line flag: 'y'
            # answers "do you have true negative controls", then blank lines
            # accept each printed default.
            try:
                chan.sendall("y\n" + "\n" * 14)
                chan.shutdown_write()
            except Exception:
                pass
            buffer = ""
            while True:
                if chan.recv_ready():
                    chunk = chan.recv(65536).decode("utf-8", "replace")
                    if not chunk:
                        break
                    buffer += chunk
                    while "\n" in buffer:
                        line, buffer = buffer.split("\n", 1)
                        sess.add(line)
                elif chan.exit_status_ready():
                    break
                else:
                    time.sleep(0.2)
            leftover = chan.recv(65536).decode("utf-8", "replace") if chan.recv_ready() else ""
            if buffer or leftover:
                sess.add(buffer + leftover)
            code = chan.recv_exit_status()
            chan.close()
            if code != 0:
                sess.add(f"\nStep failed with exit code {code}.")
                failed = True
                break
            sess.add(f"\nStep {index} finished.")
    except Exception as exc:                        # noqa: BLE001 - shown on the page
        sess.add(f"\nERROR: {exc}")
        failed = True
    sess.finish(failed)


# ------------------------------------------------------------------ rendering
STYLE = '<link rel="preconnect" href="https://fonts.googleapis.com"><link rel="preconnect" href="https://fonts.gstatic.com" crossorigin><link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Instrument+Serif:ital@0;1&family=Instrument+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap"><style>\n:root{\n  --paper:#fbf7f4; --card:#ffffff; --sunk:#f4efeb; --edge:#e4dcd6; --hair:#efe8e3;\n  --ink:#1b1517; --ink2:#5d5157; --ink3:#95868d;\n  --teal:#0f6d66; --teal-ink:#fff; --teal-wash:#e0eeec;\n  --rose:#c03a66; --rose-wash:#fbe6ed; --plum:#7a5aa6; --amber:#d1841f; --amber-wash:#fbeeda;\n  --good:#2c7c50; --warn:#96660f; --bad:#b03a42;\n  --good-bg:#e5f0e9; --warn-bg:#fbf0dc; --bad-bg:#fbe8e9;\n  --serif:"Instrument Serif",Georgia,"Times New Roman",serif;\n  --sans:"Instrument Sans",system-ui,-apple-system,"Segoe UI",Helvetica,sans-serif;\n  --mono:"IBM Plex Mono",ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;\n}\n@media (prefers-color-scheme:dark){\n  :root{\n    --paper:#15100f; --card:#1e1817; --sunk:#181312; --edge:#332a29; --hair:#261f1e;\n    --ink:#f6efec; --ink2:#b3a4a6; --ink3:#867679;\n    --teal:#62c6bc; --teal-ink:#06231f; --teal-wash:#12302c;\n    --rose:#ef7d9f; --rose-wash:#34181f; --plum:#b195d9; --amber:#e6a24e; --amber-wash:#332412;\n    --good:#5cba86; --warn:#d5a24e; --bad:#e08087;\n    --good-bg:#14281d; --warn-bg:#2d2412; --bad-bg:#2e1819;\n  }\n}\n*{box-sizing:border-box}\nhtml{-webkit-text-size-adjust:100%}\nbody{margin:0;background:var(--paper);color:var(--ink);font-family:var(--sans);\n     font-size:15.5px;line-height:1.6;-webkit-font-smoothing:antialiased}\n::selection{background:var(--rose-wash)}\n.wrap{max-width:1180px;margin:0 auto;padding:0 26px}\n\n/* ------------------------------------------------------------------- bar */\n.bar{border-bottom:1px solid var(--edge);background:var(--card)}\n.bar .wrap{display:flex;align-items:center;gap:14px;flex-wrap:wrap;padding:13px 26px}\n.logo{display:flex;align-items:center;gap:10px;font-family:var(--serif);font-size:22px;\n  letter-spacing:-.01em}\n.logo svg{display:block}\n.grow{flex:1}\n.whoami{font-family:var(--mono);font-size:11px;color:var(--ink3);text-align:right;line-height:1.4}\n.whoami b{color:var(--ink2);font-weight:500}\n.state{display:inline-flex;align-items:center;gap:7px;font-family:var(--mono);font-size:11px;\n  letter-spacing:.06em;padding:5px 11px;border:1px solid var(--edge);color:var(--ink3);\n  background:var(--sunk);border-radius:20px}\n.state i{width:6px;height:6px;border-radius:50%;background:var(--ink3)}\n.state.run{color:var(--teal);border-color:var(--teal);background:var(--teal-wash)}\n.state.run i{background:var(--teal);animation:beat 1.2s ease-in-out infinite}\n.state.ok{color:var(--good);border-color:var(--good);background:var(--good-bg)}\n.state.ok i{background:var(--good)}\n.state.bad{color:var(--bad);border-color:var(--bad);background:var(--bad-bg)}\n.state.bad i{background:var(--bad)}\n@keyframes beat{0%,100%{opacity:1}50%{opacity:.3}}\n@media (prefers-reduced-motion:reduce){.state.run i{animation:none}}\n\n/* ----------------------------------------------------------------- intro */\n.intro{border-bottom:1px solid var(--edge);\n  background:linear-gradient(180deg,var(--card) 0%,var(--paper) 100%)}\n.intro .wrap{padding:54px 26px 40px}\n.intro h1{font-family:var(--serif);font-weight:400;font-size:clamp(2.6rem,6.6vw,4.4rem);\n  line-height:1.02;letter-spacing:-.025em;margin:0 0 16px;max-width:20ch}\n.intro h1 em{font-style:italic;color:var(--rose)}\n.intro .say{margin:0;font-size:17.5px;line-height:1.6;color:var(--ink2);max-width:60ch}\n.intro .say b{color:var(--ink);font-weight:500}\n.ribbon{display:block;width:100%;max-width:720px;height:44px;margin:26px 0 0}\n\n/* the three artefacts */\n.flow{display:flex;gap:0;flex-wrap:wrap;margin:30px 0 0}\n.flow .sp{flex:1 1 210px;min-width:0;padding:14px 20px;border-top:3px solid var(--edge);\n  margin-right:2px}\n.flow .sp:nth-child(1){border-top-color:var(--teal)}\n.flow .sp:nth-child(2){border-top-color:var(--rose)}\n.flow .sp:nth-child(3){border-top-color:var(--amber)}\n.flow .sp[data-state=done]{border-top-color:var(--good)}\n.flow .sp[data-state=fail]{border-top-color:var(--bad)}\n.flow .sp[data-state=skip]{opacity:.35}\n.flow .sp[data-state=active]{background:var(--card)}\n.flow .k{font-family:var(--mono);font-size:10px;letter-spacing:.16em;text-transform:uppercase;\n  color:var(--ink3);display:block;margin-bottom:6px}\n.flow .n{font-family:var(--serif);font-size:20px;display:block;letter-spacing:-.015em;\n  margin-bottom:4px;line-height:1.15}\n.flow .f{font-family:var(--mono);font-size:10.5px;color:var(--ink3);word-break:break-all;\n  line-height:1.5;display:block}\n.flow .sp[data-state=done] .k::after{content:"  \\2713";color:var(--good)}\n\n/* ---------------------------------------------------------------- panels */\nmain{padding:40px 0 110px}\n.panel{margin:0 0 46px}\n.panel>h2{font-family:var(--serif);font-weight:400;font-size:2.05rem;letter-spacing:-.02em;\n  margin:0 0 6px;line-height:1.1;display:flex;align-items:baseline;gap:14px}\n.panel>h2 .ix{font-family:var(--mono);font-size:11px;letter-spacing:.18em;padding:4px 9px;\n  border-radius:20px;color:#fff;align-self:center}\n.panel:nth-of-type(1)>h2 .ix{background:var(--teal)}\n.panel:nth-of-type(2)>h2 .ix{background:var(--rose)}\n.panel:nth-of-type(3)>h2 .ix{background:var(--amber)}\n.panel>p.say{margin:0 0 20px;color:var(--ink2);font-size:15px;max-width:74ch}\n.sheet{background:var(--card);border:1px solid var(--edge);border-radius:14px;padding:24px 26px}\n.cols{display:grid;gap:0 28px;grid-template-columns:repeat(auto-fit,minmax(258px,1fr))}\n\n/* ----------------------------------------------------------------- forms */\nlabel{display:block;margin:0 0 20px}\n.flabel{display:block;font-size:14px;font-weight:600;color:var(--ink);letter-spacing:-.005em}\n.hint{display:block;font-size:12.5px;line-height:1.5;color:var(--ink3);margin:3px 0 7px}\n.hint code{font-family:var(--mono);font-size:11px;color:var(--rose);background:var(--rose-wash);\n  padding:1px 5px;border-radius:4px}\ninput[type=text],input[type=password],select,textarea{width:100%;padding:10px 12px;\n  border:1px solid var(--edge);border-radius:9px;background:var(--sunk);color:var(--ink);\n  font:13px/1.5 var(--mono);transition:.15s}\ninput:hover,textarea:hover,select:hover{border-color:var(--ink3)}\ninput:focus,select:focus,textarea:focus{outline:none;border-color:var(--rose);\n  background:var(--card);box-shadow:0 0 0 3px var(--rose-wash)}\ntextarea{resize:vertical;min-height:80px}\nselect{font-family:var(--sans);font-size:14px}\n.row{display:flex;gap:8px;align-items:flex-start}\n.opt{display:flex;gap:12px;align-items:flex-start;margin:0 0 16px;cursor:pointer;padding:2px 0}\n.opt input{margin:3px 0 0;accent-color:var(--rose);width:17px;height:17px;flex:none}\n.opt strong{font-size:14px;font-weight:600;display:block}\n.opt .hint{margin:2px 0 0}\n.divide{border:0;border-top:1px solid var(--hair);margin:8px 0 20px}\n\nbutton{font-family:var(--sans);cursor:pointer;transition:.15s}\n.go{background:var(--teal);color:var(--teal-ink);border:0;padding:12px 26px;border-radius:24px;\n  font-size:15px;font-weight:600;letter-spacing:-.008em}\n.go:hover{filter:brightness(1.1)}\n.go.line{background:transparent;color:var(--rose);border:1px solid var(--rose)}\n.go.line:hover{background:var(--rose-wash);filter:none}\n.go.mute{background:transparent;color:var(--ink2);border:1px solid var(--edge);font-weight:500;\n  padding:8px 16px;font-size:13px}\n.go.mute:hover{color:var(--ink);border-color:var(--ink3)}\n.pick{font-family:var(--mono);font-size:11px;padding:10px 12px;background:var(--card);\n  border:1px dashed var(--edge);color:var(--ink3);border-radius:9px;flex:none}\n.pick:hover{color:var(--rose);border-color:var(--rose);border-style:solid;background:var(--rose-wash)}\n.commit{display:flex;gap:13px;flex-wrap:wrap;align-items:center;padding:22px 0 0;\n  border-top:1px solid var(--edge)}\n.commit .hint{margin:0}\n\ndetails.more{margin:10px 0 0}\ndetails.more>summary{cursor:pointer;font-family:var(--mono);font-size:11px;letter-spacing:.1em;\n  text-transform:uppercase;color:var(--rose);list-style:none;display:inline-flex;\n  align-items:center;gap:8px;padding:6px 0}\ndetails.more>summary::-webkit-details-marker{display:none}\ndetails.more>summary::after{content:"show";font-size:10px;border-bottom:1px solid var(--rose)}\ndetails.more[open]>summary::after{content:"hide"}\ndetails.more>div{padding-top:18px;border-top:1px solid var(--hair);margin-top:6px}\n\n/* ------------------------------------------------------------------- log */\n.tty{border:1px solid var(--edge);border-radius:14px;overflow:hidden;background:var(--card)}\n.tty>header{display:flex;align-items:center;gap:11px;padding:12px 17px;background:var(--sunk);\n  border-bottom:1px solid var(--edge)}\n.tty>header .nm{font-family:var(--mono);font-size:11px;color:var(--ink2)}\n.tty>header .clock{font-family:var(--mono);font-size:11px;color:var(--ink3);margin-left:auto;\n  font-variant-numeric:tabular-nums}\npre.tap{margin:0;padding:18px;max-height:430px;overflow:auto;background:transparent;\n  font:12px/1.7 var(--mono);color:var(--ink2);white-space:pre-wrap;word-break:break-word}\n.cmd{margin:0 0 14px;padding:12px 14px;background:var(--sunk);border:1px solid var(--hair);\n  border-radius:8px;font:12px/1.6 var(--mono);color:var(--ink2);white-space:pre-wrap;\n  word-break:break-all}\npre.tap .c{color:var(--teal);font-weight:500}\npre.tap .g{color:var(--good)}\npre.tap .r{color:var(--bad)}\npre.tap .y{color:var(--warn)}\n\n/* --------------------------------------------------------------- results */\n.shelf{display:grid;gap:18px;grid-template-columns:repeat(auto-fill,minmax(228px,1fr));\n  margin-top:18px}\n.shelf a{display:block;text-decoration:none;color:inherit;border:1px solid var(--edge);\n  border-radius:12px;overflow:hidden;background:var(--card);transition:.15s}\n.shelf a:hover{border-color:var(--rose);transform:translateY(-3px)}\n.shelf .pic{aspect-ratio:16/10;background:var(--sunk);display:grid;place-items:center;\n  overflow:hidden;border-bottom:1px solid var(--hair)}\n.shelf .pic img{width:100%;height:100%;object-fit:cover;object-position:top left}\n.shelf .pic span{font-family:var(--mono);font-size:10px;letter-spacing:.12em;\n  text-transform:uppercase;color:var(--ink3)}\n.shelf .cap{padding:11px 13px;font-family:var(--mono);font-size:11px;color:var(--ink2);\n  word-break:break-all;line-height:1.45}\n\n/* --------------------------------------------------------------- notices */\n.tell{border-radius:11px;padding:13px 16px;font-size:14px;margin:0 0 18px;\n  border:1px solid var(--edge);border-left:4px solid var(--plum);background:var(--card);\n  color:var(--ink2)}\n.tell.err{background:var(--bad-bg);border-color:var(--bad);border-left-color:var(--bad);color:var(--bad)}\n.tell.warn{background:var(--warn-bg);border-color:var(--warn);border-left-color:var(--warn);color:var(--warn)}\n.tell b{font-weight:600}\nh3.cut{font-family:var(--serif);font-size:1.5rem;color:var(--ink);margin:46px 0 4px;\n  font-weight:400;letter-spacing:-.018em}\n\n/* ---------------------------------------------------------------- picker */\n.veil{position:fixed;inset:0;z-index:60;background:rgba(27,21,23,.5);display:grid;\n  place-items:center;padding:26px}\n.veil .pane{background:var(--card);border:1px solid var(--edge);border-radius:16px;\n  width:min(640px,100%);max-height:80vh;display:flex;flex-direction:column;overflow:hidden;\n  box-shadow:0 30px 70px -28px rgba(27,21,23,.5)}\n.veil>.pane>header{display:flex;align-items:center;gap:11px;padding:16px 19px;\n  border-bottom:1px solid var(--edge)}\n.veil>.pane>header strong{font-family:var(--serif);font-size:19px;font-weight:400}\n.veil .where{font-family:var(--mono);font-size:11px;color:var(--ink3);padding:10px 19px;\n  background:var(--sunk);border-bottom:1px solid var(--hair);word-break:break-all}\n.veil .rows{overflow:auto;flex:1}\n.veil>.pane>footer{padding:13px 19px;border-top:1px solid var(--edge);display:flex;gap:10px;\n  justify-content:flex-end}\n.item{display:flex;align-items:center;gap:11px;padding:10px 19px;cursor:pointer;\n  font-family:var(--mono);font-size:12px;border-bottom:1px solid var(--hair)}\n.item:hover{background:var(--rose-wash)}\n.item .g{width:14px;color:var(--ink3);flex:none;text-align:center}\n.item.d .t{color:var(--rose);font-weight:500}\n.item .t{flex:1;word-break:break-all}\n.item .s{color:var(--ink3);font-size:10.5px;font-variant-numeric:tabular-nums}\n\n/* ------------------------------------------------------------------ gate */\n.gate{min-height:100vh;display:grid;place-items:center;padding:30px;\n  background:radial-gradient(1200px 600px at 20% -10%,var(--rose-wash),transparent 60%),\n             radial-gradient(900px 500px at 100% 100%,var(--teal-wash),transparent 60%)}\n.gate .mid{width:min(452px,100%)}\n.gate .sheet{padding:34px 32px 28px}\n.gate h1{font-family:var(--serif);font-weight:400;font-size:2.3rem;letter-spacing:-.022em;\n  margin:0 0 6px;line-height:1.08}\n.gate h1 em{font-style:italic;color:var(--rose)}\n.gate .say{color:var(--ink2);font-size:14.5px;margin:0 0 24px}\n.gate .go{width:100%;margin-top:8px;padding:13px}\n.gate .fine{font-family:var(--mono);font-size:11px;color:var(--ink3);line-height:1.65;\n  margin:20px 6px 0}\n</style>'

RIBBON = '<svg class="ribbon" viewBox="0 0 720 44" role="img" aria-label="Overlapping 56-mer peptide tiles, each offset by 28 residues"><rect x="0" y="20" width="58" height="19" rx="9.5" fill="var(--teal)" opacity="0.22"/><rect x="30" y="4" width="58" height="19" rx="9.5" fill="var(--rose)" opacity="0.55"/><rect x="60" y="20" width="58" height="19" rx="9.5" fill="var(--plum)" opacity="0.22"/><rect x="90" y="4" width="58" height="19" rx="9.5" fill="var(--amber)" opacity="0.55"/><rect x="120" y="20" width="58" height="19" rx="9.5" fill="var(--teal)" opacity="0.22"/><rect x="150" y="4" width="58" height="19" rx="9.5" fill="var(--rose)" opacity="0.55"/><rect x="180" y="20" width="58" height="19" rx="9.5" fill="var(--plum)" opacity="0.22"/><rect x="210" y="4" width="58" height="19" rx="9.5" fill="var(--amber)" opacity="0.55"/><rect x="240" y="20" width="58" height="19" rx="9.5" fill="var(--teal)" opacity="0.22"/><rect x="270" y="4" width="58" height="19" rx="9.5" fill="var(--rose)" opacity="0.55"/><rect x="300" y="20" width="58" height="19" rx="9.5" fill="var(--plum)" opacity="0.22"/><rect x="330" y="4" width="58" height="19" rx="9.5" fill="var(--amber)" opacity="0.55"/><rect x="360" y="20" width="58" height="19" rx="9.5" fill="var(--teal)" opacity="0.22"/><rect x="390" y="4" width="58" height="19" rx="9.5" fill="var(--rose)" opacity="0.55"/><rect x="420" y="20" width="58" height="19" rx="9.5" fill="var(--plum)" opacity="0.22"/><rect x="450" y="4" width="58" height="19" rx="9.5" fill="var(--amber)" opacity="0.55"/><rect x="480" y="20" width="58" height="19" rx="9.5" fill="var(--teal)" opacity="0.22"/><rect x="510" y="4" width="58" height="19" rx="9.5" fill="var(--rose)" opacity="0.55"/><rect x="540" y="20" width="58" height="19" rx="9.5" fill="var(--plum)" opacity="0.22"/><rect x="570" y="4" width="58" height="19" rx="9.5" fill="var(--amber)" opacity="0.55"/><rect x="600" y="20" width="58" height="19" rx="9.5" fill="var(--teal)" opacity="0.22"/><rect x="630" y="4" width="58" height="19" rx="9.5" fill="var(--rose)" opacity="0.55"/><rect x="660" y="20" width="58" height="19" rx="9.5" fill="var(--plum)" opacity="0.22"/></svg>'

MARK = '<svg width="26" height="26" viewBox="0 0 26 26" aria-hidden="true"><rect x="1" y="8" width="15" height="10" rx="5" fill="var(--teal)" opacity=".5"/><rect x="10" y="8" width="15" height="10" rx="5" fill="var(--rose)" opacity=".9"/></svg>'


SCRIPT = '<script>\nconst $=s=>document.querySelector(s), $$=s=>[...document.querySelectorAll(s)];\n\nfunction syncPipe(){ $("#pipeBox").style.display=$("#run_pipeline").checked?"":"none"; mark(); }\nfunction syncAnn(){ const on=$("#run_annotate").checked;\n  $("#annBox").style.display=on?"":"none";\n  $("#cpmField").style.display=on?"none":"";\n  $("#seqNote").style.display=on?"":"none"; mark(); }\n$("#run_pipeline").addEventListener("change",syncPipe);\n$("#run_annotate").addEventListener("change",syncAnn);\n\n/* the three artefacts across the top, lit as the run advances */\nfunction mark(state){\n  const ann=$("#run_annotate").checked, pipe=$("#run_pipeline").checked;\n  const sp=$$("#flow .sp");\n  sp[0].dataset.state = state&&state[0] ? state[0] : (ann?"idle":"skip");\n  sp[1].dataset.state = state&&state[1] ? state[1] : "idle";\n  sp[2].dataset.state = state&&state[2] ? state[2] : (pipe?"idle":"skip");\n}\n\n/* ---- path picker ---------------------------------------------------- */\nlet pt=null, pk="file", cwd="";\ndocument.addEventListener("click",async e=>{\n  const b=e.target.closest("button.pick[data-target]"); if(!b) return;\n  pt=b.dataset.target; pk=b.dataset.kind; cwd=document.getElementById(pt).value||"";\n  $("#picker").style.display="grid"; await hop(); });\n$("#pickClose").onclick=()=>{ $("#picker").style.display="none"; };\n$("#picker").addEventListener("click",e=>{ if(e.target===$("#picker")) $("#picker").style.display="none"; });\n$("#pickUse").onclick=()=>{ document.getElementById(pt).value=cwd; $("#picker").style.display="none"; };\ndocument.addEventListener("keydown",e=>{ if(e.key==="Escape") $("#picker").style.display="none"; });\nfunction size(n){ if(n===null||n===undefined) return "";\n  const u=["B","K","M","G"]; let i=0,v=Number(n);\n  while(v>=1024&&i<3){v/=1024;i++;} return v.toFixed(v<10&&i>0?1:0)+u[i]; }\nasync function hop(){\n  const j=await (await fetch("/browse?path="+encodeURIComponent(cwd))).json();\n  if(j.error){ $("#pickList").innerHTML="<div class=\'item\'><span class=\'t\'>"+j.error+"</span></div>"; return; }\n  cwd=j.path; $("#pickPath").textContent=j.path;\n  let h="";\n  if(j.parent!==null) h+="<div class=\'item d\' data-go=\'"+j.parent+"\'><span class=\'g\'>\\u2191</span><span class=\'t\'>up one level</span></div>";\n  for(const it of j.entries){\n    if(it.dir) h+="<div class=\'item d\' data-go=\'"+j.path+"/"+it.name+"\'><span class=\'g\'>\\u25b8</span><span class=\'t\'>"+it.name+"</span></div>";\n    else if(pk==="file") h+="<div class=\'item\' data-use=\'"+j.path+"/"+it.name+"\'><span class=\'g\'>\\u00b7</span><span class=\'t\'>"+it.name+"</span><span class=\'s\'>"+size(it.size)+"</span></div>";\n  }\n  $("#pickList").innerHTML=h||"<div class=\'item\'><span class=\'t\'>nothing here</span></div>"; }\n$("#pickList").addEventListener("click",async e=>{\n  const t=e.target.closest("[data-go],[data-use]"); if(!t) return;\n  if(t.dataset.go!==undefined){ cwd=t.dataset.go; await hop(); }\n  else { document.getElementById(pt).value=t.dataset.use; $("#picker").style.display="none"; } });\n\n/* ---- log ------------------------------------------------------------ */\nlet seen=0, timer=null, t0=null, stick=true;\nconst RX=/^=== Step (\\d+)\\/(\\d+):/;\nfunction tint(line){\n  const e=document.createElement("span"); let c="";\n  if(/^=== Step/.test(line)) c="c";\n  else if(/ERROR|failed with exit code/i.test(line)) c="r";\n  else if(/Step \\d+ finished|^\\s*Wrote |identical|\\bok\\b/.test(line)) c="g";\n  else if(/^Warning:|WARNING/i.test(line)) c="y";\n  e.className=c; e.textContent=line+"\\n"; return e; }\nfunction clock(){ if(!t0) return "";\n  const s=Math.floor((Date.now()-t0)/1000);\n  return String(Math.floor(s/60)).padStart(2,"0")+":"+String(s%60).padStart(2,"0"); }\nconst tap=$("#log");\ntap.addEventListener("scroll",()=>{ stick=tap.scrollTop+tap.clientHeight>=tap.scrollHeight-24; });\nfunction say(kind,text){ const p=$("#statusPill");\n  p.className="state"+(kind?" "+kind:""); p.lastChild.textContent=text; }\nasync function poll(){\n  const j=await (await fetch("/log?since="+seen)).json();\n  if(j.lines.length){\n    if(seen===0) tap.textContent="";\n    const frag=document.createDocumentFragment();\n    for(const l of j.lines){\n      const m=RX.exec(l);\n      if(m){ const cur=parseInt(m[1],10), off=$("#run_annotate").checked?0:1;\n        const st=["idle","idle","idle"];\n        for(let i=0;i<cur-1;i++) st[i+off]="done";\n        st[cur-1+off]="active"; mark(st); }\n      frag.appendChild(tint(l));\n    }\n    tap.appendChild(frag);\n    if(stick) tap.scrollTop=tap.scrollHeight;\n  }\n  seen=j.total;\n  $("#el").textContent=clock();\n  if(j.running) say("run","working");\n  if(j.done){\n    clearInterval(timer); timer=null;\n    say(j.failed?"bad":"ok", j.failed?"stopped":"finished");\n    mark($$("#flow .sp").map(x=>x.dataset.state==="active"?(j.failed?"fail":"done"):x.dataset.state));\n    const o=await (await fetch("/outputs")).json();\n    if(o.files&&o.files.length){\n      let h="<h3 class=\'cut\'>What it made</h3><div class=\'shelf\'>";\n      for(const f of o.files){\n        const u="/file?p="+encodeURIComponent(f), nm=f.split("/").pop();\n        h+="<a href=\'"+u+"\' target=\'_blank\'><div class=\'pic\'>"\n          +(/\\.png$/i.test(f)?"<img loading=\'lazy\' src=\'"+u+"\' alt=\'\'>":"<span>"+(nm.split(".").pop()||"file")+"</span>")\n          +"</div><div class=\'cap\'>"+nm+"</div></a>";\n      }\n      $("#outputs").innerHTML=h+"</div>";\n    }\n  }\n}\nsyncAnn(); syncPipe();\nif(new URLSearchParams(location.search).get("started")){\n  seen=0; t0=Date.now(); say("run","working");\n  timer=setInterval(poll,1100); poll();\n  setInterval(()=>{ if(timer) $("#el").textContent=clock(); },1000);\n}\n</script>'


def page_setup(sess, scripts, packages, binaries, log=""):
    """Shown when the server is missing scripts or dependencies.

    Everything here runs on the person's own server over the SSH connection
    they already opened, so a colleague needs this one file and nothing else.
    Each action shows the exact command first: nothing is installed silently.
    """
    rows = []
    if scripts:
        cmd = (f"git clone {sess.repo} {sess.pipeline_dir}" if sess.repo
               else "(no repository given -- enter one at sign-in)")
        rows.append(f"""
      <div class="sheet" style="margin-bottom:16px">
        <h2>Pipeline scripts</h2>
        <p class="say">Missing from <code>{html.escape(sess.pipeline_dir)}</code>:
           {html.escape(', '.join(scripts))}</p>
        <pre class="cmd">{html.escape(cmd)}</pre>
        <form method="post" action="/setup/scripts">
          <button class="go" type="submit"{'' if sess.repo else ' disabled'}>
            Install the scripts</button>
        </form>
      </div>""")
    if packages or binaries:
        want = []
        if packages:
            want.append("Python packages: " + ", ".join(packages))
        if binaries:
            want.append("not on PATH: " + ", ".join(binaries))
        rows.append(f"""
      <div class="sheet" style="margin-bottom:16px">
        <h2>Dependencies</h2>
        <p class="say">{html.escape(' &middot; '.join(want))}</p>
        <p class="say">bootstrap.sh looks for an interpreter that already has all of
           them and only downloads something if your machine genuinely has none.
           It needs no root: everything lands in your own conda directory, about
           1&nbsp;GB. On a cluster, run this on a login node &mdash; compute nodes
           often have no outbound network.</p>
        <pre class="cmd">{html.escape(posixpath.join(sess.pipeline_dir, 'bootstrap.sh'))}</pre>
        <form method="post" action="/setup/env">
          <button class="go" type="submit">Set up the dependencies</button>
        </form>
      </div>""")
    if not rows:
        rows.append("""
      <div class="sheet" style="margin-bottom:16px">
        <h2>Ready</h2>
        <p class="say">The scripts and every dependency are present.</p>
        <form method="get" action="/"><button class="go" type="submit">Continue</button></form>
      </div>""")
    output = (f'<div class="sheet"><h2>Output</h2><pre class="cmd">{html.escape(log)}</pre></div>'
              if log else "")
    return f"""<!doctype html><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>PhIP-Seq setup</title>{STYLE}
<div class="gate"><div class="mid" style="max-width:720px">
  <div class="logo" style="margin:0 0 20px">{MARK}<span>PhIP-Seq</span></div>
  <h1 style="margin:0 0 6px">Set up {html.escape(sess.host)}</h1>
  <p class="say" style="margin:0 0 20px">Signed in as
     <code>{html.escape(sess.user)}</code>. This server needs a few things before
     the pipeline can run. Each step runs under your own account.</p>
  {''.join(rows)}
  {output}
</div></div>"""


def page_login(message="", insecure_note=""):
    err = f'<div class="tell err"><b>Could not connect.</b> {html.escape(message)}</div>' if message else ""
    warn = f'<div class="tell warn">{insecure_note}</div>' if insecure_note else ""
    return f"""<!doctype html><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>PhIP-Seq pipeline</title>{STYLE}
<div class="gate"><div class="mid">
  <div class="logo" style="margin:0 0 20px">{MARK}<span>PhIP-Seq</span></div>
  <div class="sheet">
    <h1>Sign in to the server</h1>
    <p class="say">Use your ordinary Linux account on the analysis machine.</p>
    {err}{warn}
    <form method="post" action="/login">
      <label><span class="flabel">Server</span>
        <input type="text" name="server" value="{html.escape(DEFAULT_SERVER)}"
               placeholder="analysis server hostname" required></label>
      <label><span class="flabel">Username</span>
        <input type="text" name="username" autocomplete="username" required autofocus></label>
      <label><span class="flabel">Password</span>
        <input type="password" name="password" autocomplete="current-password" required></label>
      <hr class="divide">
      <label><span class="flabel">Project folder</span>
        <span class="hint">The area you want to work in. Nothing outside it can be read, written
        or run, so pick something wide enough to reach both your data and your results.</span>
        <input type="text" name="root" value="{html.escape(DEFAULT_ROOT)}" required></label>
      <label><span class="flabel">Where the scripts live</span>
        <span class="hint">The folder holding <code>prepare_virus_matrix_v2.py</code>. Checked now,
        so a wrong path is caught here rather than half way through a run.</span>
        <input type="text" name="pipeline_dir" value="{html.escape(DEFAULT_PIPELINE_DIR)}" required></label>
      <label><span class="flabel">Python to run them with</span>
        <span class="hint">Must have pandas, numpy, matplotlib and scipy, and must find
        <code>blastp</code>, <code>tblastn</code> and <code>makeblastdb</code> on its PATH.
        A bare <code>python3</code> often does not &mdash; give the full path to the one inside
        your environment, e.g. <code>~/miniconda3/envs/phipseq/bin/python</code>. Checked now.</span>
        <input type="text" name="python_exe" value="{html.escape(DEFAULT_PYTHON)}" required></label>
      <label><span class="flabel">Scripts repository (optional)</span>
        <span class="hint">If the folder above is empty, the pipeline is installed from here
        onto your server. Leave it as it is unless you have your own fork.</span>
        <input type="text" name="repo" value="{html.escape(DEFAULT_REPO)}"
               placeholder="https://github.com/.../phipseq-pipeline.git"></label>
      <button class="go" type="submit">Sign in</button>
    </form>
  </div>
  <p class="fine">Your password opens the SSH connection and is then discarded. It is not
  stored, not written to any log, and never sent back to this page.</p>
</div></div>"""


# =============================================================================
# FIELD DESCRIPTIONS — edit these by hand.
#
#   "field_name": {
#       "label":       bold text above the input  (the only required key)
#       "help":        the grey explanatory line. HTML is allowed, so <code>x</code>
#                      and &mdash; work. Do NOT put HTML entities in "label", which
#                      is escaped — use a plain apostrophe there.
#       "value":       what the box is pre-filled with
#       "placeholder": grey example text inside an empty box
#       "pick":        "file", "dir", or omit for no browse button
#   }
#
# Leaving any key except "label" out is safe — it falls back to empty, so an edit
# cannot break the page. Restart the server to see changes.
# =============================================================================
FIELDS = {
    "plate_dirs": {
        "label": "Your phipflow results folders, one per plate",
        "help": "One path per line &mdash; the folder phipflow produced for each plate. Point at the "
            "run folder Each must contain the CPM and z-score matrices plus the peptide and "
            "sample annotation tables.",
        "placeholder": "/path/to/plate1\n/path/to/plate2",
    },
    "combined_dir": {
        "label": "Where the combined files should go",
        "placeholder": "combined",
        "pick": "dir",
    },
    "virus": {
        "label": "Short name for the virus (required)",
        "help": "Gets written into every peptide name, and is what you filter on afterwards. ",
        "placeholder": "cov2",
    },
    "segment": {
        "label": "Segment or protein name (optional)",
        "help": "Only if you are preparing one segment on its own, e.g. <code>S</code> or "
            "<code>M</code>. Leave blank when the coordinate map already names each protein.",
        "placeholder": "S",
    },
    "library": {
        "label": "All of the peptides in your library",
        "help": "The whole library&rsquo;s peptide sequences, translated to protein &mdash; one file "
            "per library (OCCUSCAN, RespScan, MICROSCAN, NOVOSCAN), in "
            "<code>Libraries/phipflow/Oligos/*_translate.fasta</code>. This is where the "
            "sequences come from: your CPM file contains numbers and names only.",
        "placeholder": "\u2026/RESPSCAN_bowtieref_phipflow_translate.fasta",
        "pick": "file",
    },
    "min_identity": {
        "label": "Minimum % identity",
        "help": "Sequence search only. Drops low-identity chance alignments. Run once "
            "with this blank, look at the spread, then choose &mdash; 96 suits "
            "SARS-CoV-2, 60 suits CCHFV M.",
        "placeholder": "60",
    },
    "ambiguity_fraction": {
        "label": "Genotype ambiguity fraction",
        "help": "References within this fraction of the top bitscore are all kept. "
            "Default 0.98.",
        "placeholder": "0.98",
    },
    "genome_dir": {
        "label": "Folder for the genome plots",
        "help": "Defaults to a <code>genome_plots</code> folder inside the run folder.",
        "placeholder": "",
    },
    "genome_length": {
        "label": "Genome or protein length",
        "help": "The x-axis maximum. 29903 for a SARS-CoV-2 genome; the protein length "
            "in amino acids for a single-protein target, e.g. 1684 for CCHFV M.",
        "placeholder": "1684",
    },
    "genome_label": {
        "label": "Reference accession",
        "help": "Printed on the x-axis, e.g. <code>NC_045512.2</code> or "
            "<code>NP_950235.1</code>.",
        "placeholder": "NP_950235.1",
    },
    "measures": {
        "label": "What to plot",
        "help": "Any of <code>pct zscore cpm q75</code>, space separated. One figure each. "
            "<code>pct</code> is breadth (how many responded), <code>zscore</code> is "
            "height (how strongly).",
        "placeholder": "pct zscore",
    },
    "z_threshold": {
        "label": "z threshold",
        "help": "Counts as positive at or above this. Default 3.5.",
        "placeholder": "3.5",
    },
    "groups": {
        "label": "Sample groups",
        "help": "One per line, as <code>Label=REGEX</code>. The regex is matched anywhere "
            "in the column name, so <code>CCHFV=_CCHFV</code> pools across buffer arms.",
        "placeholder": "Infected=^ISARIC_\nVaccinated=^DOVE_",
    },
    "search_term": {
        "label": "How to spot your virus in the peptide names",
        "help": "Text that appears in the names of your virus&rsquo;s peptides. Separate alternatives "
            "with <code>|</code> &mdash;",
        "placeholder": "EPI_ISL|severe_acute_respiratory_syndrome_coronavirus_2",
    },
    "coord_ref": {
        "label": "Reference protein(s)",
        "help": "One protein for a single-segment virus. For a virus with several proteins, put them "
            "all in one file (all ten SARS-CoV-2 ORFs, or all six Ebola NP species) and each "
            "peptide is assigned to whichever it matches best.",
        "placeholder": "\u2026/cov2_proteins.fa",
        "pick": "file",
    },
    "coord_map": {
        "label": "Turns protein positions into genome positions (optional)",
        "help": "A small tab-separated table giving each reference protein&rsquo;s start position on "
            "the genome. Needed only for whole-genome plots. Leave blank and you get amino-acid "
            "positions within the protein instead.",
        "placeholder": "\u2026/cov2_coordinate_map.tsv",
        "pick": "file",
    },
    "geno_ref": {
        "label": "Whole genomes to differentiate strains (optional)",
        "help": "Two or more complete genomes, one per strain or lineage.",
    },
    "cpm": {
        "label": "Your CPM file, straight from phipflow",
        "help": "From <code>wide_data/data_cpm.csv</code>, or an annotated / plate-combined version "
            "of it.",
        "placeholder": "cpm_combined_annotated.csv",
        "pick": "file",
    },
    "cpm_name_column": {
        "label": "Which column holds the peptide name",
        "help": "Leave blank if it is the first column. Some matrices put a plain row number first "
            "and the real name in <code>full_name</code>.",
        "placeholder": "full_name",
    },
    "cpm_drop_columns": {
        "label": "Columns that are not samples",
        "help": "Eg. Peptide ID numbers left over from Phipflow",
        "placeholder": "peptide_id",
    },
    "min_query_coverage": {
        "label": "How much of each peptide must line up",
        "help": "<code>100</code> demands the whole peptide aligns. Tiles padded with a GS linker can "
            "never reach 100 and are dropped",
        "value": "100",
    },
    "evalue": {
        "label": "BLAST significance cut-off",
        "help": "<code>10</code> is BLAST&rsquo;s own default  ",
        "value": "10",
    },
    "out_dir": {
        "label": "Where the prepared files should go",
        "help": "Created if it does not exist. The relabelled CPM matrix is written here as "
            "<code>&lt;virus&gt;[_&lt;segment&gt;]_relabelled.csv</code>.",
        "placeholder": "dove/",
        "pick": "dir",
    },
    "run_dir": {
        "label": "Where the pipeline outputs should go",
        "help": "Ranked peptides, panel tuning, the selected panel, the seropositivity report and the "
            "heatmap",
        "placeholder": "dove/run1",
        "pick": "dir",
    },
    "filter_term": {
        "label": "Which peptides the pipeline should use",
        "help": "Blank means the virus label, i.e. all of them. Add <code>AND |S|</code> to restrict "
            "the run to one protein.",
        "placeholder": "cov2",
    },
    "output_prefix": {
        "label": "Start of every output filename",
        "help": "Blank uses the virus (and segment) you entered above.",
        "placeholder": "cov2",
    },
    "mock_regex": {
        "label": "How to recognise your beads-only wells",
        "value": "mock",
        "placeholder": "control_mockIP",
    },
    "pos_regex": {
        "label": "Which samples are the positives (optional)",
        "help": "A pattern matching your known-positive sample names, e.g. <code>^DOVE_</code>.",
        "placeholder": "^DOVE_",
    },
    "neg_regex": {
        "label": "Which samples are the negatives (optional)",
        "help": "A pattern matching your known-negative or comparison samples, e.g. "
            "<code>^Healthy_control_</code>.",
        "placeholder": "^Healthy_control_",
    },
    "panel_sizes": {
        "label": "Panel sizes to try while tuning",
        "help": "The pipeline scores each size on how far real samples sit from the mocks, and "
            "suggests the best. Default at 50",
        "value": "5,10,15,20,25,30,40,50",
    },
    "cpm_thresholds": {
        "label": "CPM needed for a peptide to count as present",
        "help": "Tested during tuning. Default at 10",
        "value": "10",
    },
    "panel_choice": {
        "label": "How many peptides to actually use",
        "help": "Overrides the tuner&rsquo;s suggestion. Leave blank to accept whatever it picks.",
        "value": "50",
    },
    "cpm_choice": {
        "label": "CPM cut-off to actually use",
        "help": "Also the floor under each peptide&rsquo;s mock-derived threshold.",
        "value": "10",
    },
    "extra_args": {
        "label": "Anything else to pass to the pipeline",
        "help": "Raw flags, e.g. <code>--output-set full</code>. Must start with a dash.",
        "placeholder": "--output-set full",
    },
}


def field(name, rows=0):
    """Render one form input from its FIELDS entry. rows>0 gives a textarea."""
    spec = FIELDS[name]
    label = spec["label"]
    value = spec.get("value", "")
    placeholder = spec.get("placeholder", "")
    pick = spec.get("pick", "")
    helptext = spec.get("help", "")
    picker = (f'<button type="button" class="pick" data-target="{name}" '
              f'data-kind="{pick}">browse</button>') if pick else ""
    note = f'<span class="hint">{helptext}</span>' if helptext else ""
    if rows:
        box = (f'<textarea name="{name}" id="{name}" rows="{rows}" '
               f'placeholder="{html.escape(placeholder)}">{html.escape(value)}</textarea>')
    else:
        box = (f'<input type="text" name="{name}" id="{name}" '
               f'value="{html.escape(value)}" placeholder="{html.escape(placeholder)}">')
    return (f'<label><span class="flabel">{html.escape(label)}</span>{note}'
            f'<div class="row">{box}{picker}</div></label>')



def page_main(sess, error=""):
    err = f'<div class="tell err">{html.escape(error)}</div>' if error else ""
    modes = "".join(
        f'<option value="{m}"{" selected" if m == "mock_q75" else ""}>{m}</option>'
        for m in HEATMAP_MODES)
    return f"""<!doctype html><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>PhIP-Seq pipeline</title>{STYLE}

<div class="bar"><div class="wrap">
  <div class="logo">{MARK}<span>PhIP-Seq</span></div>
  <div class="grow"></div>
  <div class="whoami"><b>{html.escape(sess.user)}</b> on {html.escape(sess.host)}<br>
    {html.escape(sess.root)}</div>
  <span class="state" id="statusPill"><i></i><span>ready</span></span>
  <form method="post" action="/logout" style="margin:0">
    <button class="go mute" type="submit">Sign out</button></form>
</div></div>

<div class="intro">
  <div class="wrap">
    <h1>PhIP-Seq analysis in <em>three stages</em></h1>
    <p class="say">Annotate and merge phipflow results, label and genotype the outputs, then run
    a configurable pipeline that ranks and selects the most discriminating peptides and calls
    each sample.</p>
    {RIBBON}
    <div class="flow" id="flow">
      <div class="sp" data-state="idle"><span class="k">Stage one</span>
        <span class="n">Merge the plates</span>
        <span class="f">data_cpm_final_annotated.csv</span></div>
      <div class="sp" data-state="idle"><span class="k">Stage two</span>
        <span class="n">Label and genotype</span>
        <span class="f">&lt;virus&gt;_relabelled.csv</span></div>
      <div class="sp" data-state="idle"><span class="k">Stage three</span>
        <span class="n">Rank, select and call</span>
        <span class="f">selected_panel.txt &middot; seropositivity_report.csv</span></div>
    </div>
  </div>
</div>

<div class="wrap"><main>
 {err}
 <form method="post" action="/run" id="mainForm">

  <section class="panel">
   <h2><span class="ix">ONE</span>Merge the plates</h2>
   <p class="say">phipflow gives you one results folder per plate, with matrices keyed by bare
   numbers. This attaches real peptide and sample names to both the CPM and the z-score matrix,
   then merges the plates so every sample ends up in one file.</p>
   <div class="sheet">
    <label class="opt"><input type="checkbox" name="run_annotate" id="run_annotate" checked>
      <span><strong>Start from my phipflow results folders</strong>
      <span class="hint">Leave this on and the combined matrix is built here and passed straight
      to stage two. Turn it off if you already have a combined annotated CPM file.</span></span></label>
    <div id="annBox">
      <hr class="divide">
      {field("plate_dirs", rows=4)}
      {field("combined_dir")}
    </div>
   </div>
  </section>

  <section class="panel">
   <h2><span class="ix">TWO</span>Label and genotype</h2>
   <p class="say">Finds your virus in the library, works out what every peptide actually is by
   aligning its sequence, and rewrites the peptide names to carry the protein and position. The
   CPM numbers are never touched.</p>
   <p class="tell" id="seqNote">The CPM matrix comes from stage one. The translated library is
   still needed here, because that is where the peptide sequences are read from.</p>
   <div class="sheet">
    <div class="cols">
     {field("virus")}
     {field("segment")}
     {field("search_term")}
     {field("coord_ref")}
     {field("out_dir")}
    </div>
    <label class="opt"><input type="checkbox" name="sequence_search">
      <span><strong>Find the peptides by sequence instead of by name</strong>
      <span class="hint">Leave off for a library that names your virus reliably, such as
      OCCUSCAN &mdash; the search term above is then used. Turn on for one that does not,
      such as RESPSCAN: the whole library is aligned to the reference and whatever matches
      is taken as your virus. Clear the search term when you do. Aligning to a reference is
      not proof of identity, so set a minimum identity below.</span></span></label>
    <div id="cpmField" class="cols">{field("cpm")}{field("library")}</div>
    <details class="more"><summary>Coordinates, genotyping and alignment settings</summary><div>
     <div class="cols">
      {field("coord_map")}
      {field("geno_ref")}
      {field("cpm_name_column")}
      {field("cpm_drop_columns")}
      {field("min_query_coverage")}
      {field("evalue")}
      {field("min_identity")}
      {field("ambiguity_fraction")}
     </div>
     <label class="opt"><input type="checkbox" name="drop_unselected">
       <span><strong>Keep only this virus&rsquo;s peptides</strong>
       <span class="hint">Off by default: every row is kept, and peptides this run did not
       select hold on to their original names, so the matrix still works for other viruses.
       Turn on for a virus-only matrix.</span></span></label>
    </div></details>
   </div>
  </section>

  <section class="panel">
   <h2><span class="ix">THREE</span>Rank, select and call</h2>
   <p class="say">Takes the matrix just built, picks the handful of peptides that best separate
   real samples from the beads-only wells, and gives every sample a call.</p>
   <div class="sheet">
    <label class="opt"><input type="checkbox" name="run_pipeline" id="run_pipeline" checked>
      <span><strong>Run the ranking pipeline</strong>
      <span class="hint">Turn off to stop once the labelled matrix exists.</span></span></label>
    <div id="pipeBox">
     <hr class="divide">
     <div class="cols">
      {field("run_dir")}
      {field("filter_term")}
      {field("mock_regex")}
     </div>
     <details class="more"><summary>Panel tuning, sample groups and ranking mode</summary><div>
      <div class="cols">
       {field("output_prefix")}
       {field("pos_regex")}
       {field("neg_regex")}
       {field("panel_sizes")}
       {field("cpm_thresholds")}
       {field("panel_choice")}
       {field("cpm_choice")}
       <label><span class="flabel">Ranking and heatmap mode</span>
         <span class="hint">This decides which peptides enter the panel, not just how the heatmap
         is drawn. Two runs differing only here select different peptides.</span>
         <select name="heatmap_mode">{modes}</select></label>
       {field("extra_args")}
      </div>
      <label class="opt"><input type="checkbox" name="supervised">
        <span><strong>Supervised mode</strong>
        <span class="hint">Off by default. Unsupervised ranks by how far real samples sit above
        the beads-only wells and needs no labels. Supervised ranks by how well peptides separate
        your named positives from your named negatives.</span></span></label>
      <label class="opt"><input type="checkbox" name="collapse">
        <span><strong>Collapse peptides at the same position</strong>
        <span class="hint">Off by default. Keeps only the strongest where several peptides cover
        the same place. Leave it off if you want genome plots: collapsing replaces peptide names
        with bare positions, which carry no protein or coordinates.</span></span></label>
     </div></details>
    </div>
   </div>
  </section>

  <section class="panel">
   <h2><span class="ix">FOUR</span>Map the signal onto the genome</h2>
   <p class="say">Draws where along the virus the antibody signal sits, one stacked panel per
   sample group, bars coloured by protein, with the ranked panel peptides marked underneath.</p>
   <label class="opt"><input type="checkbox" name="run_genome" checked>
     <span><strong>Draw the genome plots</strong>
     <span class="hint">Uses the matrix and the panel the steps above just produced.</span></span></label>
   <div class="sheet">
    <div class="cols">
     {field("groups")}
     {field("measures")}
     {field("genome_length")}
     {field("genome_label")}
    </div>
    <details class="more"><summary>Where it lands, and the positivity threshold</summary><div>
     <div class="cols">
      {field("genome_dir")}
      {field("z_threshold")}
     </div>
     <p class="hint">Coordinates follow stage two: amino-acid positions on the reference
     protein, or genome nucleotide positions if you gave a coordinate map there.</p>
    </div></details>
   </div>
  </section>

  <div class="commit">
    <button class="go" type="submit">Run it</button>
    <button class="go line" type="submit" formaction="/preview">Show me the commands first</button>
    <span class="hint">Nothing runs until you press one of these.</span>
  </div>
 </form>

 <h3 class="cut">Progress</h3>
 <div class="tty">
  <header><span class="nm">{html.escape(sess.user)}@{html.escape(sess.host)}</span>
   <span class="clock" id="el"></span></header>
  <pre class="tap" id="log">Nothing running yet. Fill in the stages above and press Run it.</pre>
 </div>
 <div id="outputs"></div>
</main></div>

<div class="veil" id="picker" style="display:none"><div class="pane">
  <header><strong>Pick a path</strong>
    <span class="hint" style="margin:0">on {html.escape(sess.host)}</span>
    <div class="grow"></div>
    <button type="button" class="go mute" id="pickClose">Close</button></header>
  <div class="where" id="pickPath"></div>
  <div class="rows" id="pickList"></div>
  <footer><button type="button" class="go line" id="pickUse">Use this folder</button></footer>
</div></div>
{SCRIPT}"""


# ------------------------------------------------------------------- handler
class Handler(BaseHTTPRequestHandler):
    server_version = "phipseq-gui"

    def log_message(self, fmt, *args):
        if os.environ.get("PHIPSEQ_GUI_VERBOSE"):
            sys.stderr.write("%s - %s\n" % (self.address_string(), fmt % args))

    def _token(self):
        raw = self.headers.get("Cookie")
        if not raw:
            return None
        try:
            morsel = SimpleCookie(raw).get("phipseq_session")
            return morsel.value if morsel else None
        except Exception:
            return None

    def _send(self, body, status=HTTPStatus.OK, ctype="text/html; charset=utf-8", headers=()):
        data = body.encode("utf-8") if isinstance(body, str) else body
        self.send_response(status)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(data)))
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Referrer-Policy", "no-referrer")
        self.send_header("Cache-Control", "no-store")
        for key, value in headers:
            self.send_header(key, value)
        self.end_headers()
        try:
            self.wfile.write(data)
        except (BrokenPipeError, ConnectionResetError):
            pass

    def _json(self, obj, status=HTTPStatus.OK):
        self._send(json.dumps(obj), status, "application/json")

    def _form(self):
        length = int(self.headers.get("Content-Length") or 0)
        if length <= 0 or length > 1_000_000:
            return {}
        raw = self.rfile.read(length).decode("utf-8", "replace")
        return {k: v[-1] for k, v in parse_qs(raw, keep_blank_values=True).items()}

    def _redirect(self, where, headers=()):
        self.send_response(HTTPStatus.SEE_OTHER)
        self.send_header("Location", where)
        for key, value in headers:
            self.send_header(key, value)
        self.send_header("Content-Length", "0")
        self.end_headers()

    def _loopback(self):
        return self.client_address[0] in ("127.0.0.1", "::1", "localhost")

    def _insecure_note(self):
        if self._loopback():
            return ""
        if ALLOW_INSECURE:
            return ("This connection is <strong>not encrypted</strong> and you are not on "
                    "localhost. Your account password will cross the network in cleartext.")
        return ""

    def do_GET(self):
        route = urlparse(self.path)
        path, query = route.path, parse_qs(route.query)
        sess = get_session(self._token())
        if path == "/login" or sess is None:
            return self._send(page_login(insecure_note=self._insecure_note()))
        if path == "/setup":
            return self._send(page_setup(sess, *check_install(sess)))
        if path == "/":
            # Never show the run form for a server that cannot run anything --
            # the stages would fail one by one with nothing explaining why.
            scripts, packages, binaries = check_install(sess)
            if scripts or packages or binaries:
                return self._send(page_setup(sess, scripts, packages, binaries))
            return self._send(page_main(sess))
        if path == "/log":
            since = int((query.get("since") or ["0"])[0] or 0)
            lines, total, running, done, failed = sess.read(since)
            return self._json({"lines": lines, "total": total, "running": running,
                               "done": done, "failed": failed})
        if path == "/browse":
            target = (query.get("path") or [""])[0] or sess.root
            try:
                base = safe_remote(sess, target, want_dir=True)
            except FormError:
                try:
                    base = safe_remote(sess, posixpath.dirname(target) or sess.root,
                                       want_dir=True)
                except FormError as exc:
                    return self._json({"error": str(exc)})
            parent = None if base.rstrip("/") == sess.root.rstrip("/") else posixpath.dirname(base)
            try:
                return self._json({"path": base, "parent": parent,
                                   "entries": remote_listing(sess, base)})
            except IOError as exc:
                return self._json({"error": str(exc)})
        if path == "/outputs":
            if not sess.out_dir:
                return self._json({"files": []})
            pattern = " -o ".join(f"-name '*{ext}'" for ext in sorted(VIEWABLE))
            status, out = remote_run(
                sess, f"find {shlex.quote(sess.out_dir)} -maxdepth 3 -type f \\( {pattern} \\) "
                      f"| head -200")
            files = [line.strip() for line in out.splitlines() if line.strip()] if status == 0 else []
            return self._json({"files": files})
        if path == "/file":
            raw = (query.get("p") or [""])[0]
            try:
                target = safe_remote(sess, raw)
            except FormError as exc:
                return self._send(f"<p>{html.escape(str(exc))}</p>", HTTPStatus.FORBIDDEN)
            if posixpath.splitext(target)[1].lower() not in VIEWABLE:
                return self._send("<p>Not a viewable file type.</p>", HTTPStatus.FORBIDDEN)
            try:
                with sess.sftp().open(target, "rb") as handle:
                    handle.prefetch()
                    data = handle.read(20 * 1024 * 1024)
            except IOError as exc:
                return self._send(f"<p>{html.escape(str(exc))}</p>", HTTPStatus.NOT_FOUND)
            ctype = mimetypes.guess_type(target)[0] or "text/plain; charset=utf-8"
            return self._send(data, HTTPStatus.OK, ctype)
        return self._send("<p>Not found.</p>", HTTPStatus.NOT_FOUND)

    def do_POST(self):
        path = urlparse(self.path).path
        form = self._form()

        if path == "/login":
            addr = self.client_address[0]
            if not self._loopback() and not ALLOW_INSECURE:
                return self._send(page_login(
                    "Refusing to accept an account password over an unencrypted "
                    "connection from a non-local address. Use an SSH tunnel "
                    "(ssh -L <port>:127.0.0.1:<port> <this-host>), or restart the "
                    "server with --allow-insecure-credentials if you accept the risk."),
                    HTTPStatus.FORBIDDEN)
            wait = throttle(addr)
            if wait > 0:
                return self._send(page_login(
                    f"Too many failed attempts. Try again in {int(wait) + 1}s."),
                    HTTPStatus.TOO_MANY_REQUESTS)
            server = (form.get("server") or "").strip()
            username = (form.get("username") or "").strip()
            password = form.get("password") or ""
            root = (form.get("root") or DEFAULT_ROOT).strip()
            pipeline_dir = (form.get("pipeline_dir") or DEFAULT_PIPELINE_DIR).strip()
            python_exe = (form.get("python_exe") or DEFAULT_PYTHON).strip()
            port = 22
            if ":" in server:
                server, _, raw_port = server.rpartition(":")
                try:
                    port = int(raw_port)
                except ValueError:
                    return self._send(page_login("Invalid port in server address."),
                                      HTTPStatus.BAD_REQUEST)
            if not server or not username or not password:
                return self._send(page_login("Server, username and password are all required."),
                                  HTTPStatus.BAD_REQUEST)
            try:
                client = ssh_connect(server, port, username, password)
            except paramiko.AuthenticationException:
                note_failure(addr)
                return self._send(page_login("Authentication failed — check the username and password."),
                                  HTTPStatus.UNAUTHORIZED)
            except paramiko.SSHException as exc:
                note_failure(addr)
                message = str(exc)
                if "not found in known_hosts" in message.lower():
                    message = (f"Host key for {server} is not in known_hosts, so the "
                               f"connection was refused. Run `ssh {username}@{server}` "
                               f"once to record it, or restart with --auto-add-host-key.")
                return self._send(page_login(f"SSH error: {message}"), HTTPStatus.BAD_GATEWAY)
            except Exception as exc:                    # noqa: BLE001
                note_failure(addr)
                return self._send(page_login(f"Could not connect: {exc}"),
                                  HTTPStatus.BAD_GATEWAY)
            finally:
                password = None                          # drop it immediately
            clear_failures(addr)
            sess = Session(client, server, username, root.rstrip("/") or "/",
                           pipeline_dir.rstrip("/"), python_exe)
            # The root is the one thing nothing can repair -- the session is
            # confined to it, so a bad one is a dead end and we say so here.
            try:
                sess.root = sess.sftp().normalize(sess.root)
            except IOError:
                sess.close()
                return self._send(page_login(f"root directory not found: {root}"),
                                  HTTPStatus.BAD_REQUEST)
            sess.repo = (form.get("repo") or "").strip()
            # Everything else -- scripts, packages, BLAST -- is fixable from
            # here, so we admit the session and route to the setup page rather
            # than refusing the login and leaving the person nowhere.
            token = secrets.token_urlsafe(32)
            with SESSIONS_LOCK:
                SESSIONS[token] = sess
            print(f"[auth] {username}@{server} connected from {addr}", flush=True)
            scripts, packages, binaries = check_install(sess)
            landing = "/" if not (scripts or packages or binaries) else "/setup"
            return self._redirect(landing, [("Set-Cookie", f"phipseq_session={token}; Path=/; "
                                                       f"HttpOnly; SameSite=Strict")])

        if path.startswith("/setup"):
            sess = get_session(self._token())
            if sess is None:
                return self._redirect("/")
            log = ""
            if path == "/setup/scripts":
                if not sess.repo:
                    log = "No repository given. Sign out and enter one."
                else:
                    # Clone when the folder has no repo in it, pull when it does,
                    # so pressing this twice updates rather than failing.
                    target = shlex.quote(sess.pipeline_dir)
                    cmd = (f"if [ -d {target}/.git ]; then "
                           f"git -C {target} pull --ff-only; else "
                           f"mkdir -p {target} && "
                           f"git clone {shlex.quote(sess.repo)} {target}; fi")
                    print(f"[setup] {sess.user}@{sess.host}: installing scripts", flush=True)
                    status, log = remote_run(sess, cmd, timeout=300)
                    log = f"$ {cmd}\n\n{log}\n[exit {status}]"
            elif path == "/setup/env":
                boot = shlex.quote(posixpath.join(sess.pipeline_dir, "bootstrap.sh"))
                print(f"[setup] {sess.user}@{sess.host}: bootstrapping", flush=True)
                # --check first: it is read-only and usually finds an interpreter
                # that already works, which costs nothing and downloads nothing.
                status, log = remote_run(sess, f"bash {boot} --check", timeout=120)
                if status != 0:
                    status, more = remote_run(sess, f"bash {boot} </dev/null", timeout=1800)
                    log = f"{log}\n{more}"
                log = f"$ bash {boot}\n\n{log}\n[exit {status}]"
                # bootstrap prints the interpreter it settled on; adopt it so the
                # person does not have to copy it back into the form by hand.
                for line in reversed(log.splitlines()):
                    cand = line.strip()
                    if cand.endswith("/bin/python") or cand.endswith("/bin/python3"):
                        sess.python = cand
                        log += f"\n\nUsing interpreter: {cand}"
                        break
            scripts, packages, binaries = check_install(sess)
            if not (scripts or packages or binaries) and not log:
                return self._redirect("/")
            return self._send(page_setup(sess, scripts, packages, binaries, log))

        sess = get_session(self._token())
        if sess is None:
            return self._send(page_login("Session expired — please connect again."),
                              HTTPStatus.UNAUTHORIZED)
        if path == "/logout":
            drop_session(self._token())
            return self._redirect("/login",
                                  [("Set-Cookie", "phipseq_session=; Path=/; Max-Age=0")])
        if path in ("/preview", "/run"):
            try:
                steps, run_dir, matrix = build_commands(sess, form)
            except FormError as exc:
                return self._send(page_main(sess, str(exc)), HTTPStatus.BAD_REQUEST)
            if path == "/preview":
                parts = ["<h2>Commands</h2>",
                         f"<p class='muted small'>These run on "
                         f"<code>{html.escape(sess.host)}</code> as "
                         f"<code>{html.escape(sess.user)}</code>.</p>"]
                for title, cmd, workdir in steps:
                    parts.append(f"<p class='muted small'>{html.escape(title)} — in "
                                 f"<code>{html.escape(workdir)}</code></p><pre>"
                                 + html.escape(" \\\n  ".join(shlex.quote(p) for p in cmd))
                                 + "</pre>")
                parts.append(f"<p class='muted small'>Matrix passed into step 2:<br>"
                             f"<code>{html.escape(matrix)}</code></p>"
                             f"<p><a href='/'>&larr; back to the form</a></p>")
                return self._send(f"<!doctype html><meta charset='utf-8'>{STYLE}"
                                  f"<div class='wrap'>{''.join(parts)}</div>")
            if sess.running:
                return self._send(page_main(sess, "A run is already in progress."),
                                  HTTPStatus.CONFLICT)
            threading.Thread(target=execute, args=(sess, steps, matrix), daemon=True).start()
            return self._redirect("/?started=1")
        return self._send("<p>Not found.</p>", HTTPStatus.NOT_FOUND)


# ---------------------------------------------------------------------- main
def main():
    global ALLOW_INSECURE, AUTO_ADD_HOST_KEY
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--host", default="127.0.0.1",
                    help="interface to bind. Default 127.0.0.1 (this machine only)")
    ap.add_argument("--port", type=int, default=8770)
    ap.add_argument("--allow-insecure-credentials", action="store_true",
                    help="accept account passwords from non-local addresses over plain HTTP")
    ap.add_argument("--auto-add-host-key", action="store_true",
                    help="trust unknown SSH host keys on first use instead of refusing")
    ap.add_argument("--open", action="store_true")
    args = ap.parse_args()
    ALLOW_INSECURE = args.allow_insecure_credentials
    AUTO_ADD_HOST_KEY = args.auto_add_host_key

    server = ThreadingHTTPServer((args.host, args.port), Handler)
    shown = "127.0.0.1" if args.host in ("0.0.0.0", "") else args.host
    print("PhIP-seq preparation & pipeline GUI  (remote / SSH)")
    print(f"  url : http://{shown}:{args.port}/")
    print(f"  sign in with your Linux account on the analysis server")
    if args.host not in ("127.0.0.1", "localhost"):
        if ALLOW_INSECURE:
            print("\n  *** WARNING: bound to a network interface over plain HTTP with")
            print("  *** --allow-insecure-credentials. Account passwords will cross the")
            print("  *** network unencrypted. Use only on a trusted network.")
        else:
            print("\n  Bound to a network interface, but credentials are only accepted")
            print("  from localhost. Colleagues should tunnel:")
            print(f"    ssh -L {args.port}:127.0.0.1:{args.port} {os.uname().nodename}")
            print("  (or restart with --allow-insecure-credentials to accept the risk)")
    else:
        print(f"\n  remote access: ssh -L {args.port}:127.0.0.1:{args.port} {os.uname().nodename}")
    if AUTO_ADD_HOST_KEY:
        print("\n  NOTE: --auto-add-host-key is on; unknown host keys are trusted on first use.")
    print("\n  Ctrl-C to stop.\n", flush=True)
    if args.open:
        threading.Timer(0.6, lambda: webbrowser.open(f"http://{shown}:{args.port}/")).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopped.")
        with SESSIONS_LOCK:
            for sess in SESSIONS.values():
                sess.close()


if __name__ == "__main__":
    main()
