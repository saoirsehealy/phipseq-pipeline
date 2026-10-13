#!/usr/bin/env bash
# =============================================================================
# Step 0 of the PhIP-seq pipeline: get the reference files for a new virus.
#
# prepare_virus_matrix_v2.py requires --coord-ref and optionally takes
# --geno-ref, but has no way to produce either. This fetches both from NCBI by
# accession, so nothing has to be copied between servers and no reference data
# needs to live in the repository.
#
#   ./fetch_references.sh --virus cchfv --segment m \
#       --coord-acc NP_950235.1 \
#       --geno-acc  DQ133507.1,AY277672.1,DQ211638.1
#
#   ./fetch_references.sh --virus nipah --coord-acc NP_112021.1
#
# Run with no arguments and it asks for each one instead.
#
#   --virus       short name, no spaces. Used in filenames and in the labels
#                 prepare_virus_matrix_v2.py writes into peptide identifiers.
#   --segment     optional, for segmented viruses (m, s, l).
#   --coord-acc   ONE NCBI *protein* accession. This is the coordinate
#                 reference: it answers "where in the protein does this peptide
#                 sit", so every peptide gets a position on one shared axis.
#                 For a multi-ORF virus give the whole proteome instead by
#                 passing several accessions separated by commas.
#   --geno-acc    optional, comma-separated NCBI *nucleotide* accessions, one
#                 per strain or genotype you want peptides assigned to. Omit
#                 this and the run simply has no genotype column.
#
#                 Write them as LABEL=ACCESSION. The label is what appears in
#                 the genotype column of the relabelled matrix, so it should be
#                 the clade name you want to read in your results:
#
#                   --geno-acc CCHF_VI_Europe2_M=DQ211625,CCHF_I_Africa1_M=DQ211627
#
#                 This matters because NCBI does not know your clade names.
#                 DQ211625's own description is "...virus 2 strain AP92 segment
#                 M, complete sequence" -- it carries the strain but not the
#                 Europe2 assignment, which comes from the phylogenetic
#                 literature. A bare accession with no label still works, but
#                 the genotype column then reads as a strain or an accession
#                 rather than a clade.
#   --out-dir     where to write. Default: a folder named after the virus.
#
# WHAT THIS REPLACED
# ------------------
# The fetching half of setup_new_virus.sh, which also extracted and translated
# peptides and handed off to the v1 pipeline. v2 does all of that internally
# from --library, so only the NCBI half was still needed. The old script is
# kept in dormant_scripts/ for reference.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

VIRUS="" SEGMENT="" COORD_ACC="" GENO_ACCS="" OUT_DIR=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --virus)      VIRUS="${2:-}";      shift 2 ;;
    --segment)    SEGMENT="${2:-}";    shift 2 ;;
    --coord-acc)  COORD_ACC="${2:-}";  shift 2 ;;
    --geno-acc)   GENO_ACCS="${2:-}";  shift 2 ;;
    --out-dir)    OUT_DIR="${2:-}";    shift 2 ;;
    -h|--help)    sed -n '2,40p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done

# ------------------------------------------------------------------- efetch
# Resolve from PATH first so a conda env (entrez-direct, in environment.yml)
# works anywhere. The hardcoded path is zeta's install, kept only as a fallback
# -- it was the single thing stopping the old script running on another server.
if command -v efetch >/dev/null 2>&1; then
  EFETCH="$(command -v efetch)"
elif [[ -x /software/edirect-v23.4/efetch ]]; then
  EFETCH=/software/edirect-v23.4/efetch
else
  cat >&2 <<'MSG'
ERROR: efetch not found.

It comes from NCBI EDirect. The pipeline's environment includes it:
    conda env create -f environment.yml && conda activate phipseq
or install it on its own:
    conda install -c bioconda entrez-direct
MSG
  exit 1
fi

# ------------------------------------------------------------------ prompts
ask() {   # ask <varname> <prompt> <required>
  local __var="$1" __prompt="$2" __required="$3" __value=""
  [[ -n "${!__var}" ]] && return 0
  if [[ ! -t 0 ]]; then
    [[ "$__required" == yes ]] && { echo "ERROR: --${__var,,} is required" >&2; exit 2; }
    return 0
  fi
  read -rp "$__prompt" __value
  printf -v "$__var" '%s' "$__value"
  [[ "$__required" == yes && -z "${!__var}" ]] && { echo "ERROR: required" >&2; exit 2; }
  return 0
}

ask VIRUS      "Virus short name (e.g. cchfv, nipah): " yes
VIRUS="${VIRUS// /_}"
ask SEGMENT    "Segment, blank if none (m/s/l): " no
ask COORD_ACC  "NCBI PROTEIN accession for the coordinate reference: " yes
ask GENO_ACCS  "Genotyping refs as LABEL=ACCESSION, comma-separated (blank to skip): " no

SEGMENT="${SEGMENT// /}"
STEM="$VIRUS"; [[ -n "$SEGMENT" ]] && STEM="${VIRUS}_${SEGMENT}"
OUT_DIR="${OUT_DIR:-${SCRIPT_DIR}/${STEM}}"
mkdir -p "$OUT_DIR"

COORD_FA="${OUT_DIR}/${STEM}_coord_ref.fa"
GENO_FA="${OUT_DIR}/${STEM}_geno_refs.fa"

# ------------------------------------------------- coordinate reference
echo ""
echo "[1/2] Coordinate reference: $COORD_ACC"
"$EFETCH" -db protein -id "$COORD_ACC" -format fasta > "$COORD_FA"
n_coord=$(grep -c '^>' "$COORD_FA" || true)
# Count residues too, not just headers. efetch retries through transient SSL
# failures, but a download cut short still leaves a valid-looking header above
# a truncated or empty sequence, and a short coordinate reference would quietly
# misplace every peptide downstream.
aa_coord=$(grep -v '^>' "$COORD_FA" | tr -d '\n \t' | wc -c)
if [[ "$n_coord" -eq 0 || "$aa_coord" -eq 0 ]]; then
  echo "  ERROR: nothing usable downloaded ($n_coord header(s), $aa_coord residues)." >&2
  echo "         Check the accession is a PROTEIN accession, and retry." >&2
  rm -f "$COORD_FA"; exit 1
fi
echo "  -> $COORD_FA  ($n_coord sequence(s), $aa_coord residues)"
[[ "$n_coord" -gt 1 ]] && echo "     (several sequences: correct for a multi-ORF proteome, "\
"wrong if you meant a single protein)"

# ------------------------------------------------- genotyping references
if [[ -z "$GENO_ACCS" ]]; then
  echo ""
  echo "[2/2] No genotyping accessions given -- skipping. Runs will have no genotype column."
  GENO_FA=""
else
  # Split LABEL=ACC pairs into an accession list for efetch and a label map for
  # the renamer. A bare ACC with no "=" simply has no entry in the map.
  ACC_LIST="" ; LABEL_MAP="${OUT_DIR}/.${STEM}_labels.tsv" ; : > "$LABEL_MAP"
  IFS=',' read -ra _pairs <<< "$GENO_ACCS"
  for pair in "${_pairs[@]}"; do
    pair="${pair// /}"
    [[ -z "$pair" ]] && continue
    if [[ "$pair" == *=* ]]; then
      label="${pair%%=*}" ; acc="${pair#*=}"
      printf '%s\t%s\n' "${acc%%.*}" "$label" >> "$LABEL_MAP"
    else
      acc="$pair"
    fi
    ACC_LIST="${ACC_LIST:+$ACC_LIST,}$acc"
  done
  n_acc=$(echo "$ACC_LIST" | tr ',' ' ' | wc -w | tr -d ' ')
  n_lab=$(wc -l < "$LABEL_MAP" | tr -d ' ')
  echo ""
  echo "[2/2] Genotyping references: $n_acc accession(s), $n_lab labelled"
  [[ "$n_lab" -eq 0 ]] && echo "      (no LABEL=ACC pairs given -- the genotype column will" \
                               "show strains/accessions, not clade names)"
  RAW="${OUT_DIR}/${STEM}_geno_refs_raw.fa"
  "$EFETCH" -db nucleotide -id "$ACC_LIST" -format fasta > "$RAW"
  n_geno=$(grep -c '^>' "$RAW" || true)
  nt_geno=$(grep -v '^>' "$RAW" | tr -d '\n \t' | wc -c)
  echo "  downloaded: $n_geno sequence(s), $nt_geno bases"
  if [[ "$n_geno" -eq 0 || "$nt_geno" -eq 0 ]]; then
    echo "  ERROR: nothing usable downloaded. Check these are NUCLEOTIDE accessions" >&2
    rm -f "$RAW"; exit 1
  fi
  [[ "$n_geno" -ne "$n_acc" ]] && \
    echo "  WARNING: asked for $n_acc, got $n_geno -- an accession may be wrong"

  # makeblastdb rejects sequence IDs over 50 characters, and the ID is the
  # first whitespace-delimited field. Long ones are rewritten to
  # VIRUS_SEGMENT_strain_accession, which is also what ends up in the genotype
  # column of the relabelled matrix, so it has to stay readable.
  python3 - "$RAW" "$GENO_FA" "$VIRUS" "$SEGMENT" "$LABEL_MAP" <<'PY'
import re, sys

in_fa, out_fa, virus, segment, label_map = sys.argv[1:6]
MAX = 50

# accession (no version) -> the clade name the user asked for
labels = {}
for line in open(label_map):
    if line.strip():
        acc, name = line.rstrip("\n").split("\t", 1)
        labels[acc] = name

def shorten(header, idx):
    h = header.lstrip('>').strip()
    parts = h.split(None, 1)
    acc = parts[0].split('.')[0]
    desc = parts[1] if len(parts) > 1 else ""
    strain = re.search(r'strain[:\s]+(\S+)', desc, re.I)
    isolate = re.search(r'isolate[:\s]+(\S+)', desc, re.I)
    country = re.search(
        r'\b(Kenya|Uganda|China|Sudan|USA|Germany|France|Saudi|Iran|Congo|Senegal|'
        r'Nigeria|Egypt|India|Japan|Korea|Russia|Spain|Italy|Brazil|Cameroon|DRC|'
        r'Zaire|Rwanda|Gabon|Guinea|Sierra\s*Leone)\b', desc, re.I)
    tag = None
    if strain:    tag = strain.group(1)[:12].rstrip(',;')
    elif isolate: tag = isolate.group(1)[:12].rstrip(',;')
    elif country: tag = country.group(0).replace(' ', '_')[:12]
    prefix = virus.upper() + (f"_{segment}" if segment else "")
    return (f"{prefix}_{tag}_{acc}" if tag else f"{prefix}_{acc}")[:MAX]

out, used, changed, kept = [], {}, 0, 0
for line in open(in_fa):
    if not line.startswith('>'):
        out.append(line); continue
    original = line.strip()
    local_id = original.lstrip('>').split()[0]
    acc_nover = local_id.split('.')[0]
    # A user-supplied clade name always wins: it is the value that will be read
    # out of the genotype column, and NCBI has no idea what it is.
    if acc_nover in labels:
        label = f"{acc_nover}_{labels[acc_nover]}"[:MAX]
    elif len(local_id) <= MAX:
        out.append(line); kept += 1; continue
    else:
        label = shorten(original, len(used) + 1)
    base, n = label, 2
    while label in used.values():
        label = f"{base[:MAX-2]}_{n}"; n += 1
    used[original] = label
    out.append(f">{label}\n"); changed += 1

open(out_fa, 'w').writelines(out)
if changed:
    print(f"  shortened {changed} header(s) over {MAX} chars, kept {kept}")
    for orig, new in used.items():
        print(f"    {orig.lstrip('>').split()[0][:60]} -> {new}")
else:
    print(f"  all {kept} header(s) already within {MAX} chars")
PY
  rm -f "$RAW"
  echo "  -> $GENO_FA"
fi

# ----------------------------------------------------------------- next step
cat <<EOF

Done. Feed these to step two:

  python3 ${SCRIPT_DIR}/prepare_virus_matrix_v2.py \\
      --virus ${VIRUS}$([[ -n "$SEGMENT" ]] && echo " --segment ${SEGMENT}") \\
      --library  <your translated peptide library .fa> \\
      --coord-ref ${COORD_FA} \\$([[ -n "$GENO_FA" ]] && printf '\n      --geno-ref %s \\' "$GENO_FA")
      --cpm      <your CPM matrix> \\
      --output   ${OUT_DIR}/${STEM}_relabelled.csv \\
      --search-term '<regex matching your virus in the library headers>'

Use --sequence-search --min-identity <N> instead of --search-term if your
library's headers do not name viruses. The identity cutoff is virus-specific:
96 suited SARS-CoV-2, 60 suited CCHFV M.
EOF
