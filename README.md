# PhIP-seq analysis pipeline

Takes phipflow output through to called seropositives and genome-position
plots, for any peptide library and any virus.

Five stages, each usable on its own:

| | Stage | Script |
|---|---|---|
| 0 | Fetch reference files for a new virus from NCBI | `fetch_references.sh` |
| 1 | Real sample and peptide names, plates combined | `annotate_plates.py` |
| 2 | Library → virus peptides → BLAST → coordinates → genotypes | `prepare_virus_matrix_v2.py` |
| 3 | Rank, tune a panel, call, heatmap | `phipseq_rank_call_and_plot_v1.7.2.py` |
| 4 | Where the signal sits along the genome | `plot_genome.py` |

A browser front-end, `phipseq_prepare_gui.py`, chains stages 1–4 over SSH and
shows every command before it runs.

Full instructions are in **[HOW_TO_RUN.md](HOW_TO_RUN.md)**, with per-script
detail in `HOW_TO_RUN_prepare_virus_matrix_v2.md` and `HOW_TO_RUN_plot_genome.md`.

## Setup

Needs pandas, numpy, matplotlib and scipy, plus BLAST+ and NCBI EDirect on
`PATH`. `bootstrap.sh` checks all of them — the four packages, `blastp`,
`tblastn`, `makeblastdb` and `efetch` — and only downloads something if your
machine is genuinely missing one:

```bash
./bootstrap.sh            # check, and create the 'phipseq' env if needed
./bootstrap.sh --check    # report only, change nothing
```

It prints the interpreter path to use. No root required — everything installs
into your own conda directory.

## The browser GUI

`phipseq_prepare_gui.py` chains stages 1-4 over SSH: you sign in to your analysis
server with your own account, and every command is shown before it runs. It can also
install this pipeline onto a server that does not have it yet, so a colleague needs
only that one file.

```bash
pip install paramiko          # the only extra dependency, on the machine with the browser
./run_phipseq_prepare_gui.sh 8770 --open
```

Reach it at <http://127.0.0.1:8770>. It binds to loopback and refuses credentials from
anywhere else, because the browser hop is plain HTTP and would otherwise carry an
account password. From another machine, tunnel:
`ssh -L 8770:127.0.0.1:8770 <host-running-it>`.

**Set your own defaults in `site.env`**, which is not tracked, so your hostname and
paths stay local:

```bash
cp site.env.example site.env     # then edit
```

| Variable | What it is |
|---|---|
| `PHIPSEQ_SERVER` | hostname of the analysis server, e.g. `cluster.example.ac.uk`. No user, no port. Blank means type it each time. |
| `PHIPSEQ_ROOT` | the directory a session is confined to. Must contain **both** your input data and where results go — paths outside it are refused, including your peptide library. Blank means the folder these scripts sit in. |
| `PHIPSEQ_PYTHON` | interpreter on the server. `bootstrap.sh` prints the exact value to use. |
| `PHIPSEQ_REPO` | git URL the GUI installs from. Point it at your own fork if you have one. |

`phipseq_gui.py` is a second, simpler GUI covering stage 3 only, run locally.

## Bringing your own library

Nothing here is specific to one peptide library. Stage 2 makes no assumption
about tile length or overlap, and reads FASTA headers as free text, so any
translated amino-acid library works. The one choice you must make is how
peptides get selected:

- `--search-term '<regex>'` — your library's headers name the virus
- `--sequence-search --min-identity <N>` — they don't, so peptides are found by
  aligning against the coordinate reference

The identity threshold is **virus-specific and must be chosen deliberately**:
96% suited SARS-CoV-2, 60% suited CCHFV M. Inheriting the wrong one silently
loses peptides or collects unrelated ones.

## Reference data

Not in this repository. Peptide libraries are 8–90 MB each and belong in shared
storage; per-virus references are generated on demand:

```bash
./fetch_references.sh --virus cchfv --segment m \
    --coord-acc NP_950235.1 \
    --geno-acc CCHF_VI_Europe2_M=DQ211625,CCHF_I_Africa1_M=DQ211627
```

Write genotype references as `LABEL=ACCESSION`. The label becomes the value in
the genotype column of your results, and NCBI does not know your clade names —
DQ211625's own description gives the strain (AP92) but not its Europe2
assignment.

## Layout

```
bootstrap.sh          find or create a working environment
fetch_references.sh   stage 0
annotate_plates.py    stage 1, with combine_plates.py for several plates
prepare_virus_matrix_v2.py          stage 2
phipseq_rank_call_and_plot_v1.7.2.py  stage 3
plot_genome.py        stage 4, with plot_provenance.py for figure stamps
phipseq_prepare_gui.py  the browser front-end, run_*.sh to launch it
HOW_TO_RUN*.md        instructions
PROTOCOL_*.md         relabelling protocol, library-agnostic
environment.yml       every dependency, including BLAST+ and EDirect
```

This repository is the pipeline only. Experiment data, per-experiment analysis
scripts and project write-ups are deliberately excluded — `.gitignore` is an
allowlist, so shipping a new file means adding a line to it.
