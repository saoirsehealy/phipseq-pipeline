# `prepare_virus_matrix_v2.py` — how to run

Takes a translated PhIP-seq peptide library, works out which peptides belong to your
virus, maps them onto a reference protein, calls a genotype for each one, and rewrites
the peptide names in a CPM matrix so the rest of the pipeline can find them.

It replaces the chain of shell scripts in
`/datastore1/PROJECTS/PHIPSEQ/Emma/CCHFV_peptides/M/`. On CCHFV M it reproduces
Emma's output exactly — all 3,840 peptide names character-for-character identical.

For **choosing an identity threshold**, see `PROTOCOL_relabelling_any_library.md`.
This file is just how to drive the script.

```bash
SCRIPTS="/datastore1/PROJECTS/PHIPSEQ/Saoirse_practice/emma_pipeline"
LIBS="/datastore1/PROJECTS/PHIPSEQ/Libraries/phipflow/Oligos"
```

---

## What you need before you start

| Input | Flag | What it is |
|---|---|---|
| Translated library | `--library` | The whole library, not a pre-filtered extract. `$LIBS/OCCUSCAN_full_bowtieref_phipflow_translate.fasta` or the RespScan equivalent. |
| Coordinate reference | `--coord-ref` | Protein FASTA. **One** sequence for a single-protein target (CCHFV M), or the **whole proteome** for a multi-protein virus (SARS-CoV-2) — see below. |
| Genotype references | `--geno-ref` | **Multiple nucleotide** sequences. The FASTA headers become the genotype labels, so name them readably. Optional. |
| CPM matrix | `--cpm` | `data_cpm_final_annotated.csv` from the annotate step. |
| Output path | `--output` | Where the relabelled matrix goes. |

Coordinate reference must be protein, genotype references must be nucleotide. Getting
these the wrong way round is the most common mistake.

### One protein or many?

**One protein** (CCHFV M): peptides are placed on that protein and the coordinates in
the name are amino-acid positions on it. Nothing else needed.

**Many proteins** (SARS-CoV-2): put the whole proteome in `--coord-ref`. BLAST picks
the best-matching protein per peptide, so the placement also assigns the ORF. Add a
`--coord-map` to turn each protein's aa positions into genome nucleotide coordinates:

```
subject_id       label    genome_nt_start
YP_009724389.1   ORF1ab   266
YP_009724390.1   S        21563
YP_009724391.1   ORF3a    25393
...
YP_009724397.2   N        28274
```

`label` becomes the segment/ORF field in the peptide name, so with a coord map you do
not need `--segment`. `genome_nt_start` is the nucleotide position of residue 1 of that
protein. Working example: `dove/cov2_coordinate_map.tsv`.

Without a coord map the label falls back to the BLAST subject ID and coordinates stay
in amino acids. The two systems are never mixed in one file.

---

## The two ways to pick your virus's peptides

**Name-based** (`--search-term`) greps the library headers. Fast, and correct when the
library names your virus consistently — OCCUSCAN does this for CCHFV.

**Library-wide** (`--sequence-search`) BLASTs the whole library against your coordinate
reference and keeps whatever aligns. Slower, but necessary when the library names the
virus badly — RespScan files SARS-CoV-2 under bare GISAID accessions, and a name search
there misses 815 real peptides.

Use one or the other; the script refuses both together. If you don't know which your
library needs, run both and compare the resulting peptide sets. On OCCUSCAN the two
routes agreed on 99.6% of CCHFV M peptides, so either works there; on a library whose
headers are bare accessions, only the sequence search finds anything.

---

## Worked example — CCHFV M segment, OCCUSCAN

This is the verified command. It reproduces Emma's output exactly.

```bash
python3 $SCRIPTS/prepare_virus_matrix_v2.py \
  --library $LIBS/OCCUSCAN_full_bowtieref_phipflow_translate.fasta \
  --coord-ref $SCRIPTS/cchfv_m/CCHFV_M_RefSeq_NP_950235.1.fa \
  --geno-ref  $SCRIPTS/cchfv_m/CCHF_M_phylo_refs.fa \
  --cpm       $SCRIPTS/cchfv_m/data_cpm_final_annotated.csv \
  --output    cchfv_M_relabelled.csv \
  --virus cchfv --segment m \
  --search-term '(Glycoprotein|Envelopment_polyprotein).*Crimean-Congo'
```

Quote the search term — it contains `|` and `(`, which the shell will otherwise eat.

Expected output:

```
Selected 3,840 / 298,013 library peptides.
Coordinate BLAST produced 2,712 hits.
Genotype calls generated for 3,130 peptides.
  rows read        298,013
  relabelled       3,840
    of which placed 2,712
  left unchanged   294,173
```

Note the search term is not just `Crimean` — that also pulls in S-segment
nucleocapsid, giving 4,276 peptides. The protein name comes *before* `Tax=` in these
headers, which is why the regex reads protein-first.

### Same virus, library-wide instead

```bash
python3 $SCRIPTS/prepare_virus_matrix_v2.py \
  --library $LIBS/OCCUSCAN_full_bowtieref_phipflow_translate.fasta \
  --coord-ref $SCRIPTS/cchfv_m/CCHFV_M_RefSeq_NP_950235.1.fa \
  --geno-ref  $SCRIPTS/cchfv_m/CCHF_M_phylo_refs.fa \
  --cpm       $SCRIPTS/cchfv_m/data_cpm_final_annotated.csv \
  --output    cchfv_M_libwide.csv \
  --virus cchfv --segment m \
  --sequence-search --min-identity 60
```

Gives 2,673 peptides. Run it once with no `--min-identity` first to see the identity
distribution before you pick a number — see the protocol document.

---

## What the output looks like

Peptide names become:

```
<original name>|<virus>|<segment>|<genotype>|<start>-<end>
```

```
UniRef100_..._9VIRU|196-252|cchfv|m|DQ211633_CCHF_IIIa_Africa3_M|192-247
UniRef100_..._9VIRU|140-196|cchfv|m|NoGenotype|136-191
UniRef100_..._9VIRU|0-56|cchfv|m|Ambiguous6:DQ211627_...|DQ211629_...|NOCOORD
```

- `<start>-<end>` are **amino-acid positions on the coordinate reference**, unless you
  passed `--coord-map`, in which case they are genome nucleotide positions. The two are
  never mixed in one file.
- `NOCOORD` means the peptide was selected but BLAST could not place it full-length.
- `NoGenotype` means no genotype reference aligned.
- `AmbiguousN:` lists the N references within 2% of the top bitscore, pipe-separated.
  Downstream scripts expect exactly this.

Peptides your run did not select **keep their original names**, so the matrix stays
usable for other viruses. Add `--drop-unselected` if you want a virus-only matrix.

---

## All the flags

| Flag | Default | Notes |
|---|---|---|
| `--virus` | required | Goes in the name. Keep it lowercase and consistent. |
| `--segment` | none | Third field, e.g. `m`, `s`. **Pass it** — the plotting scripts expect `virus\|segment\|genotype\|coords`. |
| `--search-term` | — | Regex on the FASTA header. Mutually exclusive with `--sequence-search`. |
| `--sequence-search` | off | BLAST the whole library instead of matching names. |
| `--min-identity` | none | Only with `--sequence-search`. Drops low-identity chance alignments. |
| `--coord-map` | none | TSV of `subject_id, label, genome_nt_start` to convert aa positions to genome nucleotides. |
| `--drop-unselected` | off | Write a virus-only matrix. |
| `--evalue` | 10 | BLAST e-value. |
| `--min-query-coverage` | 100 | Peptide must align full-length. Lower it only deliberately. |
| `--ambiguity-fraction` | 0.98 | Genotype references within this fraction of the top bitscore are all kept. |

The script writes one file, the matrix named by `--output`. There is no separate log,
so paste the summary it prints into your notes.

---

## Checks before you trust the output

```bash
grep -c '^>' <coord-ref>        # must be 1
grep -c '^>' <geno-ref>         # your genotype count

head -1 <output> | cut -c1-20   # must start with: peptide,
grep -c 'NOCOORD' <output>      # unplaced; should be a minority of selected
```

Then: does the placed count look like the virus? A CCHFV M glycoprotein is ~1,684 aa,
tiled every 28 aa, so a few thousand peptides across many UniRef entries is right. A
handful means your search term or reference is wrong.

---

## Things that go wrong

**`Could not find the peptide identifier column`** — the matrix's first column is named
something unexpected. The script accepts `peptide`, `Peptide`, `peptide_id`, `ID`, `id`.
Some CPM matrices carry a capital-P `Peptide`; that works. Anything else:
`sed -i '1s/^<name>,/peptide,/' <matrix>`.

**Everything comes out `NOCOORD`** — your coordinate reference is probably the wrong
protein, or a nucleotide sequence by mistake. Check it is one protein sequence.

**Far fewer peptides than expected** — the search term. Test it before running:
```bash
grep -cE '<your regex>' $LIBS/<library>.fasta
```

**Counts differ slightly from an older run** — the C-terminal tile of each protein ends
in a trailing `X` in translated libraries (the stop codon). BLAST cannot align it, so
under `--min-query-coverage 100` those peptides used to vanish silently. The script now
strips it. This was worth 11 peptides on CCHFV M, all at the reference C-terminus.

**Shell mangles the search term** — quote it. `'(A|B).*C'`, single quotes.

---

## What comes next

The relabelled matrix goes straight into the rank/call pipeline:

```bash
python3 $SCRIPTS/phipseq_rank_call_and_plot_v1.7.2.py \
  cchfv_M_relabelled.csv \
  --filter-term cchfv \
  --call-score mock_subtracted
```

The matrix is a positional argument, not a flag.

`--filter-term` matches the virus field you set here. See `HOW_TO_RUN.md` for the rest.
