# Relabelling any virus in any library — general protocol

Uses `prepare_virus_matrix_v2.py`. Replaces the name-based extraction used up to
now. Does everything `phipseq_relabel_ccords_geno_pipeline.sh` does (coordinates by
blastp, genotype by tblastn) but blasts the **whole translated library**, so peptide
selection rests on the alignment rather than on how the library happened to name
things.

**Whole-library is the default in v2** — there is no flag to set. Supplying
`--search-term` switches back to the old name-based extraction, which is kept only
for comparison. (`--all-peptides` is accepted but does nothing; v1 required it.)

Verified against v1 on both libraries: byte-identical relabelled CSV and FASTA for
SARS-CoV-2 / RespScan and CCHFV M / OCCUSCAN.

---

## Why library-wide

Name-based extraction fails silently and the failure is invisible in the output.

Measured on RespScan / SARS-CoV-2:

- `EPI_ISL` matches 11,824 peptides, the organism string matches 859, and the two
  sets are **completely disjoint** — neither alone finds the virus
- Even using both, **815 peptides that align to the SARS-CoV-2 reference were
  missed**, because they carry SARS-CoV-1, RaTG13 or UniRef names. Conserved regions
  get filed under whichever relative the designer used.
- 55 of those 815 sit inside the SE4 epitope. The name-based matrix had SE4 tiles at
  1121 and 1149 only; library-wide recovers 1120, 1121, 1124, 1132, 1138 and 1139 —
  every tile the source publication lists.

Library-wide removes the guesswork. Cost is one BLAST of the full library — a few
minutes for OCCUSCAN (298,013 peptides) on a single core.

---

## The four inputs

| Input | Flag | Requirement |
|---|---|---|
| Translated library FASTA | `--library` | the whole thing; no pre-filtering, no `--search-term` |
| Coordinate reference | `--coord-ref` | protein FASTA; **any number of sequences** |
| Coordinate map (multi-protein only) | `--coord-map` | TSV: `subject_id, label, genome_nt_start` |
| Genotype reference | `--geno-ref` | nucleotide FASTA, **2 or more**; optional |
| CPM matrix | `--cpm` | the matrix to relabel; optional |

Libraries live in `/datastore1/PROJECTS/PHIPSEQ/Libraries/phipflow/Oligos/`:
`RESPSCAN_`, `OCCUSCAN_full_`, `NOVOSCAN_full_`, `MICROSCAN_` + `_bowtieref_phipflow_translate.fasta`.

`--annotation-table` (the run's own `data_peptide_annotation_table.csv`, oligo DNA)
is equivalent — verified identical for all 298,013 OCCUSCAN and 99,766 RespScan
peptides. Use it if you want the peptide set tied to one specific run.

---

## Step 1 — run it once with no identity filter

```bash
python3 prepare_virus_matrix_v2.py \
  --virus <name> --segment <seg> \
  --library /datastore1/PROJECTS/PHIPSEQ/Libraries/phipflow/Oligos/<LIB>_bowtieref_phipflow_translate.fasta \
  --coord-ref <reference.fa> \
  --geno-ref <lineages.fa> \
  --cpm <data_cpm_final_annotated.csv> \
  --sequence-search \
  --output <virus>_explore/<virus>_<seg>_relabelled.csv
```

`--sequence-search` is what makes this the library-wide route; without it the
script falls back to `--search-term`.

This places every peptide that aligns full-length to your reference, at any identity.
Do not use this matrix for analysis — it is the diagnostic run.

---

## Step 2 — look at the identity distribution and choose the threshold

**The threshold is virus-specific. There is no universal value.**

```bash
python3 - <<'PY'
import csv, numpy as np
v = np.array([float(r['pident']) for r in
              csv.DictReader(open('<virus>_explore/<virus>_<seg>_peptide_coords.tsv'),
                             delimiter='\t')])
print('n=%d  min %.1f  median %.1f' % (len(v), v.min(), np.median(v)))
h, e = np.histogram(v, bins=np.arange(10, 101, 5))
for c, lo in zip(h, e[:-1]):
    print('   %3d-%3d %6d %s' % (lo, lo+5, c, '#' * int(60*c/max(h.max(), 1))))
PY
```

You are looking for a **bimodal shape**: a peak of real peptides at high identity,
a long low tail of spurious cross-alignments, and a trough between them. Put the
threshold in the trough.

Two real examples, same method, completely different answers:

```
SARS-CoV-2 / RespScan (9,129 placed)      CCHFV M / OCCUSCAN (3,155 placed)
  95-100   5973  ############### peak      95-100    548  #########
  90- 95   1145  ###                       90- 95    892  ############### peak
  85- 90    635  ##                        85- 90    604  ##########
  ...                                      80- 85    273  ####
  60- 65    133  tail                      75- 80    157  ##
  20- 25      7                            65- 70     93  trough region
                                           ...
  median 98.2  ->  threshold 96            median 89.3  ->  threshold 60
```

**Set the threshold from the two groups, not from the tail alone.** For CCHFV M we
can check this directly, because the name search gives an independent answer. The
2,703 peptides both routes agree on sit at median 91% identity; the 452 the
library-wide search adds on its own sit at median 41%, and none of them is a CCHFV
entry (they are Nairobi sheep disease, Hazara, Dugbe, Tofla and other nairoviruses).
The two distributions barely overlap, so the trough is real and it is low:

```
  threshold   shared kept          extras dropped
     50%      2,665 / 2,703  98.6%    314 / 452   69%
     60%      2,640 / 2,703  97.7%    419 / 452   93%   <- best trade-off
     70%      2,537 / 2,703  93.9%    446 / 452   99%
     96%        548 / 2,703  20.3%    452 / 452  100%
```

Going above 60% buys very few extra exclusions and starts costing real peptides
fast. Read the trough off the histogram rather than carrying a number over from
another virus.

**Why they differ.** RespScan tiled SARS-CoV-2 from GISAID isolates collected up to
March 2021, all near-identical to the ancestral reference — so 96% keeps the peak.
CCHFV M lineages (Africa1-4, Asia1-2, Europe1-2) diverge substantially, and a single
lineage reference cannot be 96% identical to peptides from the others. At 96% you
keep 548 of 3,155 peptides — 17% — and lose most lineages.

**Sanity check:** a threshold that discards the tail should leave a peptide count in
the same range as a well-chosen name search. CCHFV M name-based gave 2,714;
library-wide at 60% gives 2,673 (2,640 shared + 33 extras). SARS-CoV-2 name-based at 96% gave 5,973;
library-wide at 96% gives 6,715.

---

## Step 3 — the real run

```bash
python3 prepare_virus_matrix_v2.py \
  --virus <name> --segment <seg> \
  --library <...>_translate.fasta \
  --coord-ref <reference.fa> \
  --coord-map <map.tsv> \
  --geno-ref <lineages.fa> \
  --cpm <data_cpm_final_annotated.csv> \
  --sequence-search \
  --min-identity <from step 2> \
  --output <virus>_final/<virus>_<seg>_relabelled.csv
```

`--coord-map` is optional. Without it the coordinates in the peptide names are
**amino-acid positions on the coordinate reference**, which is what Emma's script
produces; with it they are genome nucleotide positions. The two are never mixed in
one file.

By default every row of the input matrix is kept, and peptides this run did not
select keep their original names, so the matrix stays usable for other viruses.
Add `--drop-unselected` for a virus-only matrix.

---

## The references

**Coordinate reference** — the protein your positions are measured against.
One sequence for a single-protein virus; several for a multi-protein one, in which
case add `--coord-map` so amino acid positions become genome nucleotide positions:

```
subject_id       label    genome_nt_start
YP_009724390.1   S        21563
YP_009724397.2   N        28274
```

Without `--coord-map`, positions stay as amino acids and the label is `--segment`.

**Genotype reference** — one nucleotide sequence per lineage. **The FASTA header
becomes the genotype call**, so name them for the labels you want:
`>DQ211627_CCHF_I_Africa1_M`. Needs 2 or more. A peptide within 2% of the top
bitscore for several lineages is called `Ambiguous3:Africa1|Africa2|Europe1`.

Genotyping runs only on peptides that placed, so searching the whole library does
not make it slow.

---

## Output

```
<original name>|<virus>|<label>|<genotype>|<start>-<end>
```

That name is what `--filter-term` matches in the rank/call pipeline and where the
genome plots read coordinates.

The script writes one file: the relabelled CPM matrix named by `--output`. It prints
a summary of rows read, relabelled, placed and left unchanged, which is worth pasting
into your notes since there is no separate log file.

---

## Checks before trusting the output

```bash
grep -c '^>' <coord-ref>          # 1, or matches your --coord-map rows
grep    '^>' <geno-ref>           # >= 2, labelled as you want them reported
```

From the run output:

- **peptides placed** — a plausible number for the virus. Near zero means the
  reference is too distant or the wrong protein.
- **identity filter kept** — if it drops more than about half, the threshold is
  probably too strict for this virus. Go back to step 2.
- **genotyped / ambiguous** — a high ambiguous fraction means the lineage references
  are too similar to resolve 56-mers, which is information, not an error. CCHFV M:
  1,552 of 2,881 ambiguous.
- **rows relabelled** — must match what you expect in the CPM matrix.

---

## Known pitfalls

- **`makeblastdb` without `-parse_seqids`** for the coordinate database. With it,
  BLAST reports `ref|YP_009724392.1|` and no `--coord-map` row matches. The script
  already does this correctly; it matters if you run BLAST by hand.
- **Duplicate peptide names.** OCCUSCAN has 30, NOVOSCAN 1,816, disambiguated with
  `|V2`/`|V3` suffixes in the phipflow files. The script dies on a duplicate rather
  than silently keeping one.
- **`Peptide` vs `peptide`** as the CPM matrix's first column — capitalisation has
  varied between files. Use `--cpm-name-column` if the default is wrong.
- **Full query coverage is not identity.** `-qcov_hsp_perc 100` only says the
  alignment spans all 56 residues. Without `--min-identity`, peptides matching at
  ~20% across the full span are kept and placed. This is why step 2 exists.
