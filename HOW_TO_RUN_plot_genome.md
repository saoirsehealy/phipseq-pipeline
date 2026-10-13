# `plot_genome.py` — how to run

Draws where along a virus genome or protein the antibody signal sits, one stacked
panel per sample group, bars coloured by ORF.

One script in place of three. It reproduces `dove/plot_genome_maps.py` exactly (see
*Checks* below) and folds in the y-axis options of `plot_panel_genome_profiles.py`
and the mock-SD handling of `dove/plot_genome_zscore.py`. Those three are superseded
and are not part of this repository.

Input is a relabelled matrix from `prepare_virus_matrix_v2.py` — see
`HOW_TO_RUN_prepare_virus_matrix_v2.md` for how to make one.

> **Set this to your own path.** It is this install's, left here as a worked
> example; every command below uses it.
```bash
SCRIPTS="/datastore1/PROJECTS/PHIPSEQ/Saoirse_practice/emma_pipeline"
```

---

## The minimum you need

| Input | Flag | Notes |
|---|---|---|
| Relabelled matrix | `--matrix` | Full or virus-only; peptides it cannot place are dropped with a note. |
| Sample groups | `--group` | `'Label=REGEX'`, repeat for each group. **Required.** |
| Mock wells | `--mock-regex` | Default `mock`. These define the background and are excluded from groups. |

Everything else has a working default.

---

## What goes on the y axis — `--measure`

Repeatable; you get one figure per measure.

| `--measure` | y axis | Answers |
|---|---|---|
| `pct` *(default)* | % of the group positive on that peptide | **Breadth** — how many people responded |
| `zscore` | z summarised by `--stat`, with a threshold line | **Height** — how strong the response was |
| `cpm` | CPM summarised by `--stat`, with a 25–75th percentile band | Raw signal |
| `q75` | 75th percentile CPM minus the CPM floor | Matches the pipeline's `mock_q75` ranking mode |

`pct` counts a peptide positive when **z ≥ `--z-threshold` (3.5) and CPM ≥ `--cpm-floor` (10)** —
the same rule the ranking pipeline uses.

`pct` and `zscore` together are usually the most informative pair. For the DOVE data
the % axis puts vaccinated-vs-infected N reactivity at 9% vs 45%; the z axis puts it
at 0.3 vs 17.7, which is the same fact but far easier to read.

---

## The mock SD floor — read this once

z is `(CPM − mock mean) / mock SD`. Many peptides have a mock SD of **exactly zero**
— 5,314 of 9,129 in the RespScan SARS-CoV-2 set, 1,413 of 2,703 in CCHFV M — which
makes z undefined, and a near-zero SD sends it to hundreds.

So the SD is floored at the median of the non-zero mock SDs. The value used is
printed and written onto the figure.

```
mock SD floored at 2.692 (1,413 peptides at or below it, 52.3%)
```

This is more conservative than not flooring, and it only changes marginal peptides —
real epitopes are unaffected. **`--sd-floor 0` disables it**, giving a naive z and
reproducing `dove/plot_genome_maps.py` cell for cell.

---

## Choosing `--stat`

Applies to `zscore` and `cpm`.

| `--stat` | Use it |
|---|---|
| `mean` *(default)* | Comparable between groups of different size. Start here. |
| `median` | Conservative — most peptides sit near zero, because only a minority of any group responds to any one peptide. |
| `q90` | What the strongest tenth reached. Sensitive, but on a group of ~10 it is effectively the maximum. |
| `max` | **Do not use to compare groups.** A maximum grows with the number of individuals and is set by single sticky peptides. |

---

## Worked example — CCHFV M

Verified. One protein, amino-acid coordinates, no `--coord-map` upstream.

```bash
python3 $SCRIPTS/plot_genome.py \
  --matrix cchfv_M_relabelled.csv \
  --panel  cchfv_m_relabelled_selected_panel.txt \
  --out-dir cchfv_plots \
  --virus 'CCHFV M' \
  --coordinate-level aa \
  --genome-length 1684 --genome-label NP_950235.1 \
  --mock-regex mock \
  --measure pct zscore \
  --group 'CCHFV trt1=dup_trt1_CCHFV' \
  --group 'Healthy trt1=dup_trt1_neg'
```

Expected:

```
NOTE: dropping 1,137 of 3,840 peptides with no usable coordinates
  2,703 peptides | 23 mock wells | stat: mean
  mock SD floored at 2.692 (1,413 peptides at or below it, 52.3%)
  CCHFV trt1 n=8 | Healthy trt1 n=8

Peak pct per ORF:      CCHFV 25.0   Healthy 12.5
Peak zscore per ORF:   CCHFV 5.92   Healthy 1.68
```

`--genome-length` is the protein length in aa (1,684 for CCHFV M NP_950235.1).

---

## Worked example — SARS-CoV-2

Ten ORFs, genome nucleotide coordinates, built upstream with a `--coord-map`.

```bash
python3 $SCRIPTS/plot_genome.py \
  --matrix dove/prepared/cov2_relabelled.csv \
  --panel  dove/run_unsup_prepared/cov2_prep_selected_panel.txt \
  --out-dir cov2_plots \
  --virus SARS-CoV-2 \
  --coordinate-level nt \
  --genome-length 29903 --genome-label NC_045512.2 \
  --mock-regex control_mockIP \
  --measure pct zscore \
  --group 'Vaccinated=^DOVE_' \
  --group 'Infected=^ISARIC_' \
  --group 'Pre-pandemic=^Healthy_control_'
```

Note `--coordinate-level nt` and a genome length rather than a protein length. The
ORF field comes from the coord map, so all ten proteins land on one axis.

---

## What it writes

```
<prefix>_genome_<measure>.png          one per --measure, whole library
<prefix>_genome_<measure>_panel.png    same, restricted to --panel
<prefix>_genome_summary.csv            peak and median per ORF per group per measure
```

`<prefix>` defaults to `--virus` with spaces replaced. The panel figures appear only
if `--panel` is given and its peptides are found.

It also prints a peak-per-ORF table per measure, which is the table worth pasting
into your notes.

---

## Checks

The script is validated against the one it replaces:

```bash
# these two should agree cell for cell
python3 $SCRIPTS/plot_genome.py --measure pct --sd-floor 0 ...
python3 $SCRIPTS/dove/plot_genome_maps.py ...
```

On the DOVE SARS-CoV-2 set all 27 ORF × group cells match. Without `--sd-floor 0` the
large signals (S 62.5%, N 45.0%) still match; only marginal low-SD peptides differ,
which is the floor doing its job.

Sanity checks on any run:

- **Group n is what you expect.** The printed `n=` counts columns, not people. If
  samples are duplicated across buffer arms or technical replicates, n overstates
  the number of individuals.
- **The dropped count is plausible.** Unplaced (`NOCOORD`) peptides are normal;
  unannotated rows appear when you pass a full matrix rather than a virus-only one,
  which is fine.
- **Signal is clustered.** Real epitopes are blocks several peptides wide. Isolated
  single-peptide spikes scattered across the whole genome are usually background.

---

## Things that go wrong

**`ERROR: no --group matched any non-mock column`** — the regex. Groups are matched
with `re.search`, so a substring works, but check it against the real column names:
```bash
head -1 <matrix> | tr ',' '\n' | head -20
```
Mock columns are removed from groups automatically, so a pattern matching only mocks
finds nothing.

**`ERROR: no column matched --mock-regex`** — mock wells are named differently per
experiment: `trt1_mock1` in exp23, `control_mockIP_*` in DOVE.

**`ERROR: no peptide carries usable coordinates`** — the matrix was not relabelled,
or every peptide is `NOCOORD`. Check the names:
```bash
python3 -c "import csv,sys; csv.field_size_limit(10**9); r=csv.reader(open('<matrix>')); next(r); print(next(r)[0])"
```
You want `...|<virus>|<orf>|<genotype>|<start>-<end>`.

**A group spans columns with no common prefix** — `--group` takes a regex, not a
prefix, so `'CCHFV=_CCHFV'` pools `dup_trt1_CCHFV*` and `dup_trt2_CCHFV*`. Anchor
with `^` when you do want a prefix.

**Bars are invisible or fill the plot** — `--bar-width`. It defaults to 90 for nt and
1 for aa; a short protein plotted at nt width would be one solid block.

---

## Relationship to the older scripts

| Script | Status |
|---|---|
| `plot_genome.py` | Use this. |
| `dove/plot_genome_maps.py` | Superseded. Not in this repository; kept locally for older figures. |
| `plot_panel_genome_profiles.py` | Superseded. Not in this repository; kept locally for older figures. |
| `dove/plot_genome_zscore.py` | SARS-CoV-2 only — groups are hardcoded, and it crashes on `NOCOORD`. Superseded. |
