#!/usr/bin/env python3
"""Annotate each plate's phipflow output and combine the plates.

This is the step before prepare_virus_matrix.py. phipflow writes matrices whose
rows and columns are bare integers:

    data_cpm.csv          ,0,1,2,...        <- sample_id
                          0,0.0,52.6,...    <- peptide_id

Nothing downstream can use that, because the biology lives in two separate
annotation tables. This script performs the join Emma's annotate_cpm.sh does,
for the CPM *and* the z-score matrix, for every plate, then calls
combine_plates.py to merge them.

    peptide_id  -> Full_name   (from data_peptide_annotation_table.csv)
    sample_id   -> sample name (from data_sample_annotation_table.csv,
                                the fastq filename with its suffixes removed)

WHERE THE OLIGOS GO, AND WHY NOT INTO THE MATRIX
------------------------------------------------
The annotation table also carries each peptide's `oligo` (its DNA), which
setup_new_virus.sh translates to build a peptide FASTA. Those are preserved, but
in a SIDECAR file rather than as a column of the matrix.

The reason is concrete: phipseq_rank_call_and_plot treats every column that is
not named exactly `peptide` as a sample. An `oligo` column would therefore be
read as a sample whose values are DNA strings. Keeping it beside the matrix
leaves the oligos fully available, keyed by the same peptide name, without
breaking the matrix for its actual consumer.

OUTPUTS (in --out-dir)
----------------------
    <plate>_cpm_annotated.csv          per plate, peptide + real sample names
    <plate>_zscore_annotated.csv       per plate
    data_cpm_final_annotated.csv       all plates combined  <- feeds prepare_virus_matrix.py
    data_zscore_final_annotated.csv    all plates combined
    data_peptide_annotation_final.csv  peptide, oligo       <- the sidecar
    annotate_plates_log.txt            every count and check

CHECKS (all fatal, because each one silently corrupts the result)
    - the matrix columns are exactly the sample_ids in the annotation table
    - every peptide_id in the matrix has an annotation row
    - Full_name is unique, or the later relabelling join is ambiguous
    - sample names are unique within and across plates, or columns collide
    - all plates share the same peptide_id -> Full_name mapping, i.e. the same
      library; combining plates from different libraries is meaningless

Usage:
    python3 annotate_plates.py plate1/results plate2/results --out-dir combined/
    python3 annotate_plates.py plate1 plate2 plate3 --out-dir combined/
    python3 annotate_plates.py <one plate> --out-dir out/          # nothing to combine
"""
import argparse
import csv
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
COMBINE = HERE / "combine_plates.py"
FASTQ_SUFFIXES = ("_merged_fil.fastq", "_merged.fastq", ".fastq", ".fq.gz", ".fastq.gz")

csv.field_size_limit(10_000_000)


def die(message):
    raise SystemExit(f"ERROR: {message}")


def sample_name(row):
    """Readable sample name from the fastq path, as annotate_cpm.sh does."""
    raw = (row.get("fastq_filepath") or "").split("/")[-1].strip()
    for suffix in FASTQ_SUFFIXES:
        if raw.endswith(suffix):
            raw = raw[: -len(suffix)]
            break
    return raw or f"sample_{row.get('sample_id', '?')}"


def read_peptides(path):
    peptides, oligos = {}, {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        # Matched case-insensitively: phipflow writes peptide_id/Full_name/oligo,
        # but the capitalisation has varied between versions.
        found = {str(c).strip().lower(): c for c in (reader.fieldnames or [])}
        columns = {}
        for column in ("peptide_id", "full_name"):
            if column not in found:
                die(f"{path.name} has no '{column}' column (found {reader.fieldnames})")
            columns[column] = found[column]
        oligo_col = found.get("oligo")
        for row in reader:
            pid = row[columns["peptide_id"]].strip()
            peptides[pid] = row[columns["full_name"]].strip()
            if oligo_col:
                oligos[pid] = (row.get(oligo_col) or "").strip()
    if not peptides:
        die(f"{path.name} contains no peptides")
    seen = {}
    for pid, name in peptides.items():
        if name in seen:
            die(f"Full_name is not unique in {path.name}: {name!r} used by "
                f"peptide_id {seen[name]} and {pid}. The relabelling join would be "
                f"ambiguous, so this must be fixed upstream.")
        seen[name] = pid
    return peptides, oligos, bool(oligo_col)


def read_samples(path):
    names, statuses = {}, {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        if "sample_id" not in (reader.fieldnames or []):
            die(f"{path.name} has no 'sample_id' column")
        for row in reader:
            sid = row["sample_id"].strip()
            names[sid] = sample_name(row)
            statuses[sid] = (row.get("control_status") or "").strip()
    if not names:
        die(f"{path.name} contains no samples")
    flipped = {}
    for sid, name in names.items():
        if name in flipped:
            die(f"two samples in {path.name} resolve to the same name {name!r} "
                f"(sample_id {flipped[name]} and {sid})")
        flipped[name] = sid
    return names, statuses


def annotate_matrix(matrix_path, peptides, sample_names, out_path, log):
    """Stream the matrix, replacing integer ids with names. Low memory."""
    with open(matrix_path, newline="") as fin, open(out_path, "w", newline="") as fout:
        reader = csv.reader(fin)
        writer = csv.writer(fout)
        header = next(reader, None)
        if not header:
            die(f"{matrix_path.name} is empty")
        ids = [h.strip() for h in header[1:]]
        missing = [i for i in ids if i not in sample_names]
        if missing:
            die(f"{matrix_path.name} has {len(missing)} column(s) with no annotation "
                f"row, first few: {missing[:5]}")
        if set(ids) != set(sample_names):
            extra = sorted(set(sample_names) - set(ids))
            die(f"{matrix_path.name} covers {len(ids)} samples but the annotation "
                f"table lists {len(sample_names)}; not in the matrix: {extra[:5]}")
        # Per-plate files carry peptide_id + full_name, which is the layout
        # combine_plates.py --annotated-cpm requires. The COMBINED file is then
        # rewritten peptide-first, matching the existing data_cpm_final_annotated.csv
        # convention and prepare_virus_matrix.py's default of "name in column 1".
        writer.writerow(["peptide_id", "full_name"] + [sample_names[i] for i in ids])
        rows = unmapped = 0
        for record in reader:
            if not record:
                continue
            rows += 1
            pid = record[0].strip()
            name = peptides.get(pid)
            if name is None:
                unmapped += 1
                continue
            writer.writerow([pid, name] + record[1:])
        if unmapped:
            die(f"{matrix_path.name}: {unmapped} of {rows} rows had a peptide_id with "
                f"no annotation row. Refusing to write a partial matrix.")
    log(f"    {out_path.name}: {rows:,} peptides x {len(ids)} samples")
    return rows, len(ids)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", metavar="PLATE",
                    help="each plate's phipflow output: the run folder, its results/ "
                         "folder, or wide_data/ itself — whichever you have")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--keep-per-plate", action="store_true",
                    help="keep the per-plate annotated files (default: keep them anyway, "
                         "this flag is accepted for clarity)")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    log_lines = []

    def log(text=""):
        print(text)
        log_lines.append(text)

    plates = []
    for raw in args.inputs:
        given = Path(raw).expanduser().resolve()
        if not given.is_dir():
            die(f"not a directory: {given}")
        # Accept whatever the user points at: the phipflow run directory, its
        # results/ folder, or wide_data/ itself. A normal experiment gives one
        # results folder per plate, so that is the natural thing to hand over.
        wide = None
        for candidate in (given, given / "wide_data", given / "results" / "wide_data"):
            if (candidate / "data_cpm.csv").is_file():
                wide = candidate
                break
        if wide is None:
            die(f"no data_cpm.csv found in {given}, {given}/wide_data or "
                f"{given}/results/wide_data — is this a phipflow output folder?")
        needed = {"cpm": wide / "data_cpm.csv",
                  "zscore": wide / "data_zscore.csv",
                  "pep": wide / "data_peptide_annotation_table.csv",
                  "samp": wide / "data_sample_annotation_table.csv"}
        for key, path in needed.items():
            if not path.is_file():
                die(f"{wide} is missing {path.name}")
        plates.append((given, wide, needed))

    log(f"Annotating {len(plates)} plate(s)")
    reference = None
    cpm_outputs, z_outputs = [], []
    oligo_table, has_oligo_any = {}, False

    for index, (given, wide, files) in enumerate(plates, start=1):
        label = given.name if given.name not in ("wide_data", "results") else (
            wide.parent.parent.name if wide.name == "wide_data" else wide.name)
        label = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in label) or f"plate{index}"
        log(f"\n  [{index}] {given}")
        if wide != given:
            try:
                log(f"      using {wide.relative_to(given)}/")
            except ValueError:
                log(f"      using {wide}")
        peptides, oligos, has_oligo = read_peptides(files["pep"])
        sample_names, statuses = read_samples(files["samp"])
        beads = sum(1 for v in statuses.values() if v == "beads_only")
        log(f"    {len(peptides):,} peptides | {len(sample_names)} samples "
            f"({beads} beads_only) | oligo column: {'yes' if has_oligo else 'NO'}")
        if reference is None:
            reference = peptides
        elif peptides != reference:
            differing = sum(1 for k in set(peptides) | set(reference)
                            if peptides.get(k) != reference.get(k))
            die(f"plate {index} has a different peptide_id -> Full_name mapping "
                f"({differing:,} entries differ). All plates must use the same library.")
        if has_oligo:
            has_oligo_any = True
            for pid, seq in oligos.items():
                oligo_table[peptides[pid]] = seq

        cpm_out = out_dir / f"{label}_cpm_annotated.csv"
        z_out = out_dir / f"{label}_zscore_annotated.csv"
        annotate_matrix(files["cpm"], peptides, sample_names, cpm_out, log)
        annotate_matrix(files["zscore"], peptides, sample_names, z_out, log)
        cpm_outputs.append(cpm_out)
        z_outputs.append(z_out)

    # sample names must not collide between plates
    seen = {}
    for path in cpm_outputs:
        with open(path, newline="") as handle:
            for name in next(csv.reader(handle))[2:]:
                if name in seen:
                    die(f"sample name {name!r} appears in both {seen[name]} and "
                        f"{path.name}; combining would overwrite one of them")
                seen[name] = path.name

    final_cpm = out_dir / "data_cpm_final_annotated.csv"
    final_z = out_dir / "data_zscore_final_annotated.csv"

    def to_peptide_first(src, dest):
        """Rewrite peptide_id+full_name -> a single 'peptide' column."""
        with open(src, newline="") as fin, open(dest, "w", newline="") as fout:
            reader, writer = csv.reader(fin), csv.writer(fout)
            header = next(reader)
            try:
                i_name = header.index("full_name")
            except ValueError:
                die(f"{Path(src).name} has no full_name column to promote")
            drop = {header.index("peptide_id"), i_name} if "peptide_id" in header else {i_name}
            keep = [i for i in range(len(header)) if i not in drop]
            writer.writerow(["peptide"] + [header[i] for i in keep])
            n = 0
            for row in reader:
                if row:
                    writer.writerow([row[i_name]] + [row[i] for i in keep])
                    n += 1
        return n

    if len(plates) == 1:
        log(f"\n  single plate: no combining needed")
        for src, dest, what in ((cpm_outputs[0], final_cpm, "CPM"),
                                (z_outputs[0], final_z, "z-score")):
            n = to_peptide_first(src, dest)
            log(f"    {what}: wrote {dest.name} ({n:,} peptides)")
    else:
        log(f"\n  Combining {len(plates)} plates with combine_plates.py")
        for sources, target, what in ((cpm_outputs, final_cpm, "CPM"),
                                      (z_outputs, final_z, "z-score")):
            staged = out_dir / f"_combined_{what.replace('-', '')}_with_ids.csv"
            cmd = [sys.executable, str(COMBINE), "--annotated-cpm",
                   "--output-file", str(staged)] + [str(p) for p in sources]
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                sys.stderr.write(result.stdout + result.stderr)
                die(f"combine_plates.py failed on the {what} matrices")
            n = to_peptide_first(staged, target)
            staged.unlink()
            log(f"    {what}: wrote {target.name} ({n:,} peptides)")

    # sidecar: the oligos, keyed by the same peptide name as the matrix
    sidecar = out_dir / "data_peptide_annotation_final.csv"
    if has_oligo_any:
        with open(sidecar, "w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["peptide", "oligo"])
            for name, seq in oligo_table.items():
                writer.writerow([name, seq])
        log(f"    oligos: wrote {sidecar.name} ({len(oligo_table):,} peptides)")
    else:
        log("    oligos: no 'oligo' column in any annotation table; sidecar not written")

    with open(final_cpm, newline="") as handle:
        head = next(csv.reader(handle))
    log(f"\nFinal matrix: {final_cpm}")
    log(f"  column 1 = {head[0]!r} | {len(head) - 1} sample columns")
    log(f"\nNext: prepare_virus_matrix.py --cpm {final_cpm}")
    (out_dir / "annotate_plates_log.txt").write_text("\n".join(log_lines) + "\n")


if __name__ == "__main__":
    main()
