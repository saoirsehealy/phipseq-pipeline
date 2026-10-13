#!/usr/bin/env python3
"""Combine PhIPflow ``wide_data`` outputs from two or more plates.

All input directories must use the same peptide library. Numeric sample IDs
are rebased so the combined annotation table and numeric matrices agree.

Example:
    python combine_plates.py --output-dir /analysis/combined/wide_data \
        /analysis/plate1/wide_data /analysis/plate2/wide_data
"""

import argparse
from pathlib import Path

import pandas as pd


NUMERIC_FILES = (
    "data_counts.csv",
    "data_cpm.csv",
    "data_zscore.csv",
    "data_enrichment.csv",
)
ANNOTATED_FILES = (
    "data_counts_annotated.tsv",
    "data_cpm_annotated.tsv",
    "data_zscore_annotated.tsv",
)
SAMPLE_ANNOTATION = "data_sample_annotation_table.csv"
PEPTIDE_ANNOTATION = "data_peptide_annotation_table.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Combine PhIPflow wide_data outputs from two or more plates."
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        type=Path,
        metavar="INPUT",
        help=(
            "Each plate's PhIPflow wide_data directory, or (with --annotated-cpm) "
            "each annotated CPM CSV; input order determines sample-column order."
        ),
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        help="Directory in which to write the combined wide_data files.",
    )
    parser.add_argument(
        "--annotated-cpm",
        action="store_true",
        help="Merge annotated CPM CSVs directly instead of complete wide_data directories.",
    )
    parser.add_argument(
        "--output-file",
        type=Path,
        help="Output CSV path for --annotated-cpm mode.",
    )
    return parser.parse_args()


def require_inputs(plate_dirs: list[Path]) -> None:
    if len(plate_dirs) < 2:
        raise ValueError("Provide at least two plate directories to combine.")

    required = (SAMPLE_ANNOTATION, PEPTIDE_ANNOTATION, *NUMERIC_FILES)
    for directory in plate_dirs:
        if not directory.is_dir():
            raise ValueError(f"Plate directory does not exist: {directory}")
        missing = [filename for filename in required if not (directory / filename).is_file()]
        if missing:
            raise ValueError(
                f"{directory} is missing required file(s): {', '.join(missing)}"
            )


def files_present_for_every_plate(
    filenames: tuple[str, ...], plate_dirs: list[Path]
) -> tuple[str, ...]:
    """Return optional matrix types available for all plates, reporting skips."""
    available = []
    for filename in filenames:
        present = [directory for directory in plate_dirs if (directory / filename).is_file()]
        if len(present) == len(plate_dirs):
            available.append(filename)
        elif present:
            absent = [str(directory) for directory in plate_dirs if directory not in present]
            print(f"Skipping {filename}: absent from {', '.join(absent)}")
    return tuple(available)


def validate_same_peptides(plate_dirs: list[Path]) -> pd.DataFrame:
    reference = pd.read_csv(plate_dirs[0] / PEPTIDE_ANNOTATION)
    for directory in plate_dirs[1:]:
        candidate = pd.read_csv(directory / PEPTIDE_ANNOTATION)
        if not reference.equals(candidate):
            raise ValueError(
                "Peptide annotation tables differ; the plates cannot be safely combined: "
                f"{plate_dirs[0]} and {directory}"
            )
    return reference


def load_samples(plate_dirs: list[Path]) -> tuple[pd.DataFrame, list[dict[str, str]]]:
    combined = []
    id_maps = []
    next_id = 0

    for directory in plate_dirs:
        samples = pd.read_csv(directory / SAMPLE_ANNOTATION)
        if "sample_id" not in samples.columns:
            raise ValueError(f"{directory / SAMPLE_ANNOTATION} has no sample_id column")
        if samples["sample_id"].duplicated().any():
            raise ValueError(f"{directory / SAMPLE_ANNOTATION} contains duplicate sample IDs")

        old_ids = samples["sample_id"].astype(str)
        new_ids = [str(sample_id) for sample_id in range(next_id, next_id + len(samples))]
        id_maps.append(dict(zip(old_ids, new_ids)))
        samples = samples.copy()
        samples["sample_id"] = [int(sample_id) for sample_id in new_ids]
        combined.append(samples)
        next_id += len(samples)

    return pd.concat(combined, ignore_index=True), id_maps


def validate_and_relabel_numeric_columns(
    matrix: pd.DataFrame, id_map: dict[str, str], source: Path
) -> pd.DataFrame:
    columns = [str(column) for column in matrix.columns]
    if len(columns) != len(set(columns)):
        raise ValueError(f"{source} contains duplicate numeric sample columns")
    expected = set(id_map)
    found = set(columns)
    if found != expected:
        missing = sorted(expected - found)
        unexpected = sorted(found - expected)
        details = []
        if missing:
            details.append(f"missing sample IDs: {', '.join(missing)}")
        if unexpected:
            details.append(f"unexpected columns: {', '.join(unexpected)}")
        raise ValueError(f"{source} does not match its sample annotation table ({'; '.join(details)})")
    matrix.columns = [id_map[column] for column in columns]
    return matrix


def combine_numeric_file(
    filename: str, plate_dirs: list[Path], id_maps: list[dict[str, str]], output_dir: Path
) -> None:
    matrices = []
    reference_index = None
    for directory, id_map in zip(plate_dirs, id_maps):
        source = directory / filename
        matrix = pd.read_csv(source, index_col=0)
        matrix = validate_and_relabel_numeric_columns(matrix, id_map, source)
        if reference_index is None:
            reference_index = matrix.index
        elif not reference_index.equals(matrix.index):
            raise ValueError(f"Peptide rows differ in {filename}: {plate_dirs[0]} and {directory}")
        matrices.append(matrix)

    combined = pd.concat(matrices, axis=1)
    combined.to_csv(output_dir / filename)
    print(f"  {filename}: {combined.shape[0]} peptides x {combined.shape[1]} samples")


def combine_annotated_file(filename: str, plate_dirs: list[Path], output_dir: Path) -> None:
    matrices = []
    reference_index = None
    seen_samples: set[str] = set()
    for directory in plate_dirs:
        matrix = pd.read_csv(directory / filename, sep="\t", index_col=0)
        if reference_index is None:
            reference_index = matrix.index
        elif not reference_index.equals(matrix.index):
            raise ValueError(f"Peptide rows differ in {filename}: {plate_dirs[0]} and {directory}")
        duplicates = seen_samples.intersection(matrix.columns)
        if duplicates:
            raise ValueError(
                f"Duplicate annotated sample name(s) in {filename}: {', '.join(sorted(duplicates))}"
            )
        seen_samples.update(matrix.columns)
        matrices.append(matrix)

    combined = pd.concat(matrices, axis=1)
    combined.to_csv(output_dir / filename, sep="\t")
    print(f"  {filename}: {combined.shape[0]} peptides x {combined.shape[1]} samples")


def combine_annotated_cpm_files(cpm_files: list[Path], output_file: Path) -> None:
    """Merge already-annotated CPM CSVs while preserving peptide metadata."""
    if len(cpm_files) < 2:
        raise ValueError("Provide at least two annotated CPM files to combine.")
    if any(not filename.is_file() for filename in cpm_files):
        missing = [str(filename) for filename in cpm_files if not filename.is_file()]
        raise ValueError(f"Annotated CPM file(s) not found: {', '.join(missing)}")

    frames = [pd.read_csv(filename) for filename in cpm_files]
    required = {"peptide_id", "full_name"}
    for filename, frame in zip(cpm_files, frames):
        missing = required - set(frame.columns)
        if missing:
            raise ValueError(f"{filename} is missing required column(s): {', '.join(sorted(missing))}")

    reference = frames[0]
    reference_ids = reference["peptide_id"]
    reference_names = reference["full_name"]
    sample_matrices = []
    seen_samples: set[str] = set()
    for filename, frame in zip(cpm_files, frames):
        if not reference_ids.equals(frame["peptide_id"]) or not reference_names.equals(frame["full_name"]):
            raise ValueError(
                f"Peptide IDs or full names differ in {filename}; cannot safely combine the matrices"
            )
        sample_columns = [
            column for column in frame.columns if column not in {"peptide_id", "full_name"}
        ]
        duplicates = seen_samples.intersection(sample_columns)
        if duplicates:
            raise ValueError(
                f"Duplicate sample name(s) across CPM files: {', '.join(sorted(duplicates))}"
            )
        seen_samples.update(sample_columns)
        sample_matrices.append(frame[sample_columns])

    combined = pd.concat(sample_matrices, axis=1)
    combined.insert(0, "peptide_id", reference_ids)
    combined["full_name"] = reference_names
    output_file.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output_file, index=False)
    print(f"Combined {len(cpm_files)} annotated CPM files: "
          f"{combined.shape[0]} peptides x {len(seen_samples)} samples")
    print(f"Wrote: {output_file}")


def main() -> None:
    args = parse_args()
    if args.annotated_cpm:
        if args.output_dir is not None:
            raise ValueError("Use --output-file, not --output-dir, with --annotated-cpm")
        if args.output_file is None:
            raise ValueError("--annotated-cpm requires --output-file")
        combine_annotated_cpm_files(
            [filename.resolve() for filename in args.inputs], args.output_file.resolve()
        )
        return

    if args.output_dir is None:
        raise ValueError("--output-dir is required unless using --annotated-cpm")
    if args.output_file is not None:
        raise ValueError("--output-file can only be used with --annotated-cpm")
    plate_dirs = [directory.resolve() for directory in args.inputs]
    output_dir = args.output_dir.resolve()
    require_inputs(plate_dirs)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Combining {len(plate_dirs)} plates:")
    for directory in plate_dirs:
        print(f"  - {directory}")

    peptide_annotation = validate_same_peptides(plate_dirs)
    peptide_annotation.to_csv(output_dir / PEPTIDE_ANNOTATION, index=False)

    samples, id_maps = load_samples(plate_dirs)
    samples.to_csv(output_dir / SAMPLE_ANNOTATION, index=False)
    print(f"Sample annotation: {len(samples)} samples")
    print(f"Peptide annotation: {len(peptide_annotation)} peptides")

    for filename in NUMERIC_FILES:
        combine_numeric_file(filename, plate_dirs, id_maps, output_dir)
    for filename in files_present_for_every_plate(ANNOTATED_FILES, plate_dirs):
        combine_annotated_file(filename, plate_dirs, output_dir)

    print(f"\nDone. Combined files written to: {output_dir}")


if __name__ == "__main__":
    try:
        main()
    except ValueError as error:
        raise SystemExit(f"Error: {error}")
