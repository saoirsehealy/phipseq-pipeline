# PhIP-seq pipeline — how to run

The whole process, phipflow output to genome plot.

```bash
SCRIPTS="/datastore1/PROJECTS/PHIPSEQ/Saoirse_practice/emma_pipeline"
LIBS="/datastore1/PROJECTS/PHIPSEQ/Libraries/phipflow/Oligos"
```

| Step | Script | What it does |
|---|---|---|
| 1 | phipflow | one run per plate → `wide_data/` |
| 2 | `annotate_plates.py` | real sample and peptide names; merges plates |
| 3 | `prepare_virus_matrix_v2.py` | picks your virus's peptides, maps coordinates, genotypes, relabels |
| 4 | `run_phipseq_rank_call_and_plot_v1.7.2.sh` | ranks peptides, tunes a panel, calls samples, heatmap |
| 5 | `plot_genome.py` | where the signal sits along the genome |

Steps 2 and 3 each have their own detailed guide:
`HOW_TO_RUN_prepare_virus_matrix_v2.md`, `HOW_TO_RUN_plot_genome.md`.

---

## Step 0 — Reference files for a new virus

Only needed the first time you work on a virus. Skip it if the coord-ref and
geno-ref files already exist.

Step 3 requires `--coord-ref` and optionally takes `--geno-ref`, but cannot
produce either. `fetch_references.sh` downloads both from NCBI by accession, so
no reference data has to be copied between machines.

```bash
./fetch_references.sh \
    --virus cchfv --segment m \
    --coord-acc NP_950235.1 \
    --geno-acc CCHF_VI_Europe2_M=DQ211625,CCHF_I_Africa1_M=DQ211627,\
CCHF_IVb_Asia2_M=DQ211629,CCHF_V_Europe1_M=DQ211631,\
CCHF_IVa_Asia1_M=DQ211632,CCHF_IIIa_Africa3_M=DQ211633,\
CCHF_II_Africa2_M=DQ211637,CCHF_IIIb_Africa4_M=KX013451
```

Run it with no arguments and it asks for each value instead.

**Write the genotype refs as `LABEL=ACCESSION`.** The label becomes the value in
the genotype column of the relabelled matrix, so it should be the clade name you
want to read in your results. NCBI does not know your clade names — DQ211625's
own description is *"...virus 2 strain AP92 segment M, complete sequence"*, which
carries the strain but not the Europe2 assignment. Give a bare accession and the
genotype column shows a strain or an accession instead of a clade.

`--coord-acc` is a **protein** accession, `--geno-acc` are **nucleotide**
accessions. The command above regenerates `cchfv_m/CCHF_M_phylo_refs.fa`
byte-for-byte, all eight sequences.

> An older `setup_new_virus.sh` did this too, plus peptide extraction for the
> previous route. Step 3 now does the peptide half itself, so that script is
> retired and `fetch_references.sh` is the only way in.

---

## Step 1 — phipflow

Run it separately on each plate; the mock correction is per-plate. Each run gives a
`wide_data/` folder whose matrices are indexed by bare integers:

```
data_cpm.csv                      ,0,1,2,...      <- sample_id
data_zscore.csv                   0,0.0,52.6,...  <- peptide_id
data_sample_annotation_table.csv  sample_id  -> sample name
data_peptide_annotation_table.csv peptide_id -> Full_name
```

Nothing downstream can use those integers, which is what Step 2 is for.

---

## Step 2 — Annotate and combine the plates

```bash
python3 $SCRIPTS/annotate_plates.py plate1/results plate2/results --out-dir combined/

# one plate only
python3 $SCRIPTS/annotate_plates.py <plate> --out-dir out/
```

Each argument can be the run folder, its `results/`, or `wide_data/` itself —
whichever you have. It joins `peptide_id → Full_name` and `sample_id → sample name`
for the **CPM and the z-score** matrix, for every plate, then calls
`combine_plates.py` to merge them. You do not run `combine_plates.py` yourself.

Outputs:

```
data_cpm_final_annotated.csv       <- the one Step 3 takes
data_zscore_final_annotated.csv
```

Both have `peptide` as the first column, which is what Step 3 expects.

Two things it does deliberately. The `oligo` sequences go to a **sidecar file**, not
a column, because the rank/call script treats every column not named exactly
`peptide` as a sample — an `oligo` column would become a sample full of DNA strings.
And it **fails** if `Full_name` is not unique or the plates disagree on the
`peptide_id → Full_name` mapping, because the Step 3 join would be ambiguous.

---

## Step 3 — Relabelling and genotyping

`prepare_virus_matrix_v2.py` turns a bare peptide ID into an informative name:

```
pep_00123   ->   pep_00123|CCHFV|M|Uganda|12345-12510
```

so peptides can be grouped by protein, placed on the genome, and read in the context
of viral diversity.

### The five inputs

| Input | Flag | Notes |
|---|---|---|
| Translated library FASTA | `--library` | the whole library, e.g. OCCUSCAN or RESPSCAN |
| Viral **protein** reference | `--coord-ref` | one protein, or the whole proteome for a multi-ORF virus |
| Coordinate map *(optional)* | `--coord-map` | turns protein aa positions into genome nucleotide positions |
| Genotype **nucleotide** references *(optional)* | `--geno-ref` | headers become the genotype labels |
| CPM matrix | `--cpm` | `data_cpm_final_annotated.csv` from Step 2 |

**`--coord-ref` ≠ `--geno-ref`.** Coord-ref is protein and answers *where does this
peptide map*. Geno-ref is nucleotide and answers *which genotype is it compatible
with*. Genotyping is optional; without `--geno-ref` every peptide is `NoGenotype`.

### Two ways to pick the peptides

**Well-labelled library (OCCUSCAN) — name search, the default.** Matches the FASTA
headers with a **case-insensitive regular expression** (not literal grep), builds a
temporary virus FASTA, and BLASTs only those peptides.

```bash
python3 $SCRIPTS/prepare_virus_matrix_v2.py \
  --library $LIBS/OCCUSCAN_full_bowtieref_phipflow_translate.fasta \
  --search-term '(Glycoprotein|Envelopment_polyprotein).*Crimean-Congo' \
  --coord-ref cchfv_m/CCHFV_M_RefSeq_NP_950235.1.fa \
  --geno-ref  cchfv_m/CCHF_M_phylo_refs.fa \
  --cpm       data_cpm_final_annotated.csv \
  --output    cchfv_M_relabelled.csv \
  --virus cchfv --segment m
```

Quote the search term — it contains `|` and `(`. Be specific: plain `Crimean` also
matches S-segment nucleocapsid (4,276 peptides instead of 3,840).

**Poorly-labelled library (RESPSCAN) — sequence search.** BLASTs the entire
translated library against the coordinate reference and keeps what aligns.

```bash
python3 $SCRIPTS/prepare_virus_matrix_v2.py \
  --library $LIBS/RESPSCAN_bowtieref_phipflow_translate.fasta \
  --sequence-search --min-identity 96 \
  --coord-ref dove/work/cov2_proteins.fa \
  --coord-map dove/cov2_coordinate_map.tsv \
  --cpm       data_cpm_final_annotated.csv \
  --output    cov2_relabelled.csv \
  --virus cov2
```

Do **not** read `--sequence-search` as "these are definitely SARS-CoV-2 peptides".
It means they align to the SARS-CoV-2 reference. Run once with no `--min-identity`,
look at the identity distribution, then pick the threshold — 96% suits SARS-CoV-2,
60% suits CCHFV M. See `PROTOCOL_relabelling_any_library.md`.

### Genotype calling

Highest-scoring reference is the candidate; references within **98% of the top
bitscore** are retained. One qualifying reference gives that reference, several give
`AmbiguousN:ref1|ref2|...`, none gives `NoGenotype`.

### Output

One file, named by `--output`. **Use a `.csv` extension** — a `.tsv` extension writes
a tab-separated file and Step 4 asks for a CSV. Non-selected peptides keep their
original names so the matrix stays usable for other viruses; add `--drop-unselected`
for a virus-only matrix (Step 5's genome plots are happy either way).

---

## Step 4 — Ranking, panel tuning and calling

```bash
cd <output directory>
$SCRIPTS/run_phipseq_rank_call_and_plot_v1.7.2.sh
```

The wrapper calls `phipseq_rank_call_and_plot_v1.7.2.py`. Feed it the relabelled
matrix from Step 3.

### What the wrapper asks

| Question | What it changes |
|---|---|
| Path to input CPM matrix CSV | the relabelled matrix from Step 3 |
| Peptide filter term | which peptides enter. `cov2` takes them all; `cov2 AND \|S\|` restricts to one protein |
| Output prefix | the start of every output filename |
| Output directory | where the run lands |
| **Heatmap / ranking mode** | the most consequential answer — see below |
| ELISA CSV | optional overlay; Enter to skip |
| Use location collapse? | keep only the strongest peptide per position. **Default N here.** Turn it off if you want genome plots |
| Additional arguments | raw flags passed through; leave empty |
| Compare with a second run? | runs `compare_runs_coverage.R` |

Then the Python script asks:

| Question | What it changes |
|---|---|
| True negative controls? | `y` if you have real negatives. `n` switches to common-virus mode, thresholds from the mocks alone |
| Mock sample regex | how it finds the beads-only wells; `mock` works for both libraries here |
| Supervised mode? | what the ranking optimises. Start unsupervised |
| Positive / negative regex | define the comparison in supervised mode; labels only in unsupervised |
| Panel sizes to test | default `5,10,15,20,25,30,40,50` |
| CPM cut-offs to test | how much signal counts as present; 10 in all runs here |
| Minimum hits for a positive call | default 5 |
| Include mocks in the heatmap | display only |

### The heatmap / ranking mode

`--heatmap-mode` sets **both** the heatmap transform **and** the ranking method, so
it decides which peptides enter the panel and therefore every downstream figure.

For each peptide, `v` is the sample's CPM, `f` the CPM floor, and `μ`, `σ`, `Q75`,
`max`, `median` come from that peptide's mock wells.

| # | Name | Signal kept | The question it asks |
|---|---|---|---|
| 1 | `raw_cpm` | `v` if `v ≥ f` | detection above the floor only; no background model |
| 2 | `zscore` | `(v − μ)/σ` if `v ≥ f` | distance above the mock mean in mock SDs. Scored 0 when `σ = 0` |
| 3 | `minus_mean` | `max(0, v − μ − f)` | excess over the expected background |
| 4 | `positive_mock_median` | `max(0, v − B − f)` | excess over the median of the mock wells above zero (`B`; the mean if none are) |
| 5 | `mock_q75` | `max(0, v − Q75 − f)` | excess over the 75th percentile of the background |
| 6 | `mock_q100` | `max(0, v − max − f)` | excess over the highest mock value for that peptide |
| 7 | `mock_q100_plus_median` | `max(0, v − max − median − f)` | excess over the highest mock value plus a median mock value |
| 8 | `mean_3_5sd` | `v` if `v > max(f, μ + 3.5σ)` | parametric outlier test against the mocks |

**Mode 7 is the Python script's internal default**, but **the wrapper's prompt
defaults to 5 / `mock_q75`** — pressing Enter gives 5, not 7. Type the mode you want.

### How peptides are scored

*Unsupervised* — enrichment above the beads-only wells; it never sees your labels,
so it finds the most enriched peptides, **not the most discriminating**. A peptide
many sera bind non-specifically ranks highly without meaning anything.

```
breadth     = number of non-mock samples with signal, + 1
fold_change = (mean non-mock + 1) / (mean mock + 1)
margin      = max(mean non-mock − mean mock, 0) + 1
score       = breadth × fold_change × margin
```

*Supervised* — scores each peptide on how well it separates the samples you named as
positive from those you named as negative, minus a penalty for peptides the mocks
bind.

### Panel tuning

Eight sizes are tried and scored on how far the samples sit from the mocks, in hits
and total signal, divided by the spread among mocks. The script defaults to its own
tuned suggestion; override with `--panel-choice`.

**50 is the largest size tested**, so a tuned answer of 50 means the grid ran out,
not that 50 is optimal. Widen with the "panel sizes to test" answer if you want more.

### How a sample gets called

Per peptide, from the mocks: `threshold = max(CPM floor, mock mean + 3.5 × mock SD)`.
Then per sample:

```
hits < 3   OR  score < t_low      -> Non-reactive
hits >= 5  AND score >= t_high    -> Reactive
anything else                     -> Indeterminate
```

with `t_low = median + 3σ` and `t_high = median + 5σ` of the negatives.

> **The calling step is currently broken.** Under the default score the thresholds
> sit above the maximum achievable value, so runs come back with **zero Reactive**
> even for confirmed infections. `--call-score mock_subtracted` works around it.
> Treat `Non-reactive` from a default run as "not yet called", not as a result.

### In a browser instead

`phipseq_gui.py` runs the same pipeline from a web page, every prompt on one form.

```bash
$SCRIPTS/run_phipseq_gui.sh        # then open http://127.0.0.1:8765
ssh -L 8765:127.0.0.1:8765 user@this-host      # if working remotely
```

"Preview command" shows the command line before you commit. Leaving **CPM cut-off**
or **number of top peptides** blank accepts the pipeline's tuned suggestion, which is
only worked out part-way through the run.

### Outputs

```
<prefix>_ranked_peptides.csv          every peptide, scored
<prefix>_panel_tuning.csv             the size/cut-off grid
<prefix>_selected_panel.txt           the chosen panel  -> Step 5
<prefix>_seropositivity_report.csv    one row per sample, call and hit count
<prefix>_run_settings.csv             what you actually answered
<prefix>_heatmap_<mode>.png
```

---

## Step 5 — Genome plots

```bash
python3 $SCRIPTS/plot_genome.py \
  --matrix cov2_relabelled.csv \
  --panel  cov2_selected_panel.txt \
  --out-dir genome_plots \
  --virus SARS-CoV-2 \
  --coordinate-level nt \
  --genome-length 29903 --genome-label NC_045512.2 \
  --mock-regex control_mockIP \
  --measure pct zscore \
  --group 'Vaccinated=^DOVE_' \
  --group 'Infected=^ISARIC_' \
  --group 'Pre-pandemic=^Healthy_control_'
```

For a single-protein target use `--coordinate-level aa` and the protein length
(1,684 for CCHFV M).

`--measure` picks the y axis: `pct` (% of the group positive — breadth), `zscore`
(how strong — height), `cpm`, `q75`. Give several and you get one figure per
measure. `--group` takes a **regex**, so `'CCHFV=_CCHFV'` pools across buffer arms.

**The panel file must come from a run with location collapse OFF.** Collapsed panels
hold bare position labels (`24977-25144`) that carry no ORF or coordinates, and the
plot cannot use them.

Full detail in `HOW_TO_RUN_plot_genome.md`.

### Other plots

```bash
python3 $SCRIPTS/plot_phipseq_seroprevalence.py --exp-dir <run dir> --viruses cchfv
```

`plot_genome.py` replaces the three earlier plotting scripts
(`plot_panel_genome_profiles.py`, `dove/plot_genome_maps.py`,
`dove/plot_genome_zscore.py`). Those remain in the working directory for
reproducing older figures but are not part of this repository.

---

## Things to watch out for

**The filter term must match the labels Step 3 added.** `Africa2 AND |s|` returns
nothing if the genotype or segment labels are not in the peptide names.

**Filter terms match as substrings.** `--filter-term cchfv` also catches peptides
whose UniRef ID happens to contain `CCHFV`, such as S-segment nucleocapsid entries
named `RepID=NCAP_CCHFV`. Use `'|cchfv|'` with the pipes, or build the matrix with
`--drop-unselected`.

**The wrapper's ranking-mode default is 5, not 7**, even though the menu above it
labels 7 as the pipeline default.


**Sample group labels** in the seroprevalence scripts come from patterns in the
sample names (the `label_sample()` function). New sample types fall into `Other` and
go uncoloured until you add them there and to `GROUP_COLOURS`.
