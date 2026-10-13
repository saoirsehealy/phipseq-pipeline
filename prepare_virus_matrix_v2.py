#!/usr/bin/env python3
# =============================================================================
# PhIP-seq virus preparation + coordinate/genotype relabelling
# =============================================================================
#
# PURPOSE
# -------
# Prepare a virus-specific CPM matrix for:
#
#     phipseq_rank_call_and_plot_v1.7.2.py
#
# PEPTIDE SELECTION
# -----------------
# Default:
#     --search-term
#
#     Select peptides from the translated library using their FASTA headers.
#
# Optional:
#     --sequence-search
#
#     BLAST the entire translated library against the target-virus protein
#     reference and use matching peptides as the candidate virus peptide set.
#
#     This is intended for poorly named/inconsistently annotated libraries
#     such as cases where the library headers do not reliably identify the
#     virus.
#
# DOWNSTREAM ANNOTATION
# ---------------------
# Selected peptides are:
#
#     1. BLASTP'd against the coordinate protein reference
#     2. Mapped to genome nucleotide coordinates if --coord-map is supplied
#     3. Optionally genotyped using TBLASTN
#
# GENOTYPING
# ----------
# Follows Emma-style 98% top-bitscore logic:
#
#     top bitscore = best genotype/reference hit
#     retain every reference >= 98% of top bitscore
#
#     one reference:
#         genotype = reference name
#
#     multiple references:
#         genotype = AmbiguousN:ref1|ref2|...
#
# FINAL PEPTIDE NAME
# ------------------
#
#     <original>|<virus>|<protein>|<genotype>|<start>-<end>
#
# Example:
#
#     pep_00123|CCHFV|M|Uganda|12345-12510
#
# Ambiguous example:
#
#     pep_00124|CCHFV|M|Ambiguous2:Uganda|Nigeria|12500-12665
#
# IMPORTANT
# ---------
# The input library should be the TRANSLATED peptide FASTA containing the
# actual peptide sequences used by the PhIP-seq experiment.
#
# =============================================================================

import argparse
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd


# =============================================================================
# FASTA FUNCTIONS
# =============================================================================

def read_fasta(path):
    """
    Read FASTA.

    Returns a list of dictionaries containing:

        header
        id
        sequence

    The ID is the first whitespace-delimited field in the FASTA header.

    Example:

        >pep_00001 some annotation here

    becomes:

        id = pep_00001
        header = pep_00001 some annotation here
    """

    records = []

    header = None
    sequence = []

    with open(path) as fh:

        for line in fh:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if header is not None:

                    full_sequence = "".join(sequence)

                    records.append({
                        "header": header,
                        "id": header.split()[0],
                        "sequence": full_sequence,
                    })

                header = line[1:]
                sequence = []

            else:

                sequence.append(line)

        if header is not None:

            full_sequence = "".join(sequence)

            records.append({
                "header": header,
                "id": header.split()[0],
                "sequence": full_sequence,
            })

    return records


def write_fasta(records, path):
    """
    Write FASTA records.

    Records can be dictionaries containing:
        id
        sequence

    or:
        header
        sequence
    """

    with open(path, "w") as fh:

        for record in records:

            fasta_id = record.get(
                "id",
                record["header"].split()[0]
            )

            # Translated libraries render the stop codon of a C-terminal tile as a
            # trailing X (or *). BLAST cannot align it, so with -qcov_hsp_perc 100
            # the whole peptide is dropped. Strip it so the last tile of each
            # protein is placed like every other tile.
            sequence = record["sequence"].rstrip("X*x")

            fh.write(f">{fasta_id}\n")
            fh.write(f"{sequence}\n")


# =============================================================================
# NAME-BASED PEPTIDE SELECTION
# =============================================================================

def select_by_name(records, search_term):
    """
    Select peptides whose FASTA header matches search_term.

    search_term is treated as a regular expression.

    Multiple alternatives can be supplied with |.

    Example:

        EPI_ISL|severe_acute_respiratory_syndrome_coronavirus_2
    """

    pattern = re.compile(
        search_term,
        re.IGNORECASE
    )

    selected = []

    for record in records:

        if pattern.search(record["header"]):

            selected.append(record)

    return selected


# =============================================================================
# COMMAND RUNNER
# =============================================================================

def run_command(command):
    """
    Run a command and stop if it fails.
    """

    print("\nRunning:")
    print(" ".join(str(x) for x in command))

    subprocess.run(
        command,
        check=True
    )


# =============================================================================
# BLASTP
# =============================================================================

def run_blastp(
    query_fasta,
    reference_fasta,
    output_tsv,
    evalue=10,
    min_query_coverage=100,
    max_target_seqs=1,
):
    """
    BLAST peptide queries against a protein reference.

    Output columns:

        qseqid
        sseqid
        sstart
        send
        pident
        length
        qlen
        qcovs
        bitscore
    """

    command = [
        "blastp",

        "-query",
        str(query_fasta),

        "-subject",
        str(reference_fasta),

        "-out",
        str(output_tsv),

        "-evalue",
        str(evalue),

        "-max_target_seqs",
        str(max_target_seqs),

        "-max_hsps",
        "1",

        "-qcov_hsp_perc",
        str(min_query_coverage),

        "-outfmt",
        (
            "6 "
            "qseqid "
            "sseqid "
            "sstart "
            "send "
            "pident "
            "length "
            "qlen "
            "qcovs "
            "bitscore"
        ),
    ]

    run_command(command)


# =============================================================================
# SEQUENCE-BASED LIBRARY DISCOVERY
# =============================================================================

def sequence_search_library(
    library_fasta,
    reference_fasta,
    output_fasta,
    evalue=10,
    min_query_coverage=100,
    min_identity=None,
):
    """
    BLAST the entire translated library against the target-virus protein
    reference.

    This is NOT the default selection method.

    It is intended for libraries where peptide names/headers are incomplete
    or unreliable.

    The first/highest-scoring target is retained for each peptide because
    this is a discovery step, not the final coordinate/genotype assignment.
    """

    with tempfile.TemporaryDirectory() as tmpdir:

        blast_output = Path(tmpdir) / "library_discovery.tsv"

        run_blastp(
            query_fasta=library_fasta,
            reference_fasta=reference_fasta,
            output_tsv=blast_output,
            evalue=evalue,
            min_query_coverage=min_query_coverage,
            max_target_seqs=1,
        )

        if (
            not blast_output.exists()
            or blast_output.stat().st_size == 0
        ):

            print(
                "\nNo peptides matched the target reference "
                "during sequence-based discovery."
            )

            write_fasta(
                [],
                output_fasta
            )

            return []

        hits = pd.read_csv(
            blast_output,
            sep="\t",
            header=None,
            names=[
                "qseqid",
                "sseqid",
                "sstart",
                "send",
                "pident",
                "length",
                "qlen",
                "qcovs",
                "bitscore",
            ],
        )

    # Optional identity threshold
    if min_identity is not None:

        hits = hits[
            hits["pident"] >= min_identity
        ]

    hit_ids = set(
        hits["qseqid"].astype(str)
    )

    library_records = read_fasta(
        library_fasta
    )

    selected = [
        record
        for record in library_records
        if record["id"] in hit_ids
    ]

    write_fasta(
        selected,
        output_fasta
    )

    print(
        f"\nSequence-based discovery selected "
        f"{len(selected):,} / "
        f"{len(library_records):,} library peptides."
    )

    return selected


# =============================================================================
# COORDINATE MAP
# =============================================================================

def load_coord_map(path):
    """
    Load coordinate map.

    Required columns:

        subject_id
        label
        genome_nt_start
    """

    coord_map = pd.read_csv(
        path,
        sep="\t"
    )

    required = {
        "subject_id",
        "label",
        "genome_nt_start",
    }

    missing = required - set(
        coord_map.columns
    )

    if missing:

        raise ValueError(
            "Coordinate map is missing required columns: "
            + ", ".join(sorted(missing))
        )

    return coord_map


# =============================================================================
# COORDINATE MAPPING
# =============================================================================

def map_coordinates(
    blast_df,
    coord_map,
):
    """
    Convert protein-relative BLAST positions to genome nucleotide coordinates.

    BLAST coordinates are 1-based amino-acid coordinates.

    genome_nt_start is assumed to be the nucleotide position corresponding
    to amino-acid position 1 of the reference protein.
    """

    coord_lookup = coord_map.set_index(
        "subject_id"
    )

    results = []

    for _, row in blast_df.iterrows():

        query_id = str(
            row["qseqid"]
        )

        subject_id = str(
            row["sseqid"]
        )

        start_aa = min(
            int(row["sstart"]),
            int(row["send"])
        )

        end_aa = max(
            int(row["sstart"]),
            int(row["send"])
        )

        if subject_id not in coord_lookup.index:

            results.append({
                "qseqid": query_id,
                "sseqid": subject_id,
                "label": "NOCOORD",
                "start_aa": start_aa,
                "end_aa": end_aa,
                "genome_start": pd.NA,
                "genome_end": pd.NA,
                "pident": row["pident"],
                "bitscore": row["bitscore"],
            })

            continue

        reference = coord_lookup.loc[
            subject_id
        ]

        # Protect against duplicate subject IDs in the map.
        if isinstance(
            reference,
            pd.DataFrame
        ):

            reference = reference.iloc[0]

        genome_nt_start = int(
            reference["genome_nt_start"]
        )

        genome_start = (
            genome_nt_start
            + (start_aa - 1) * 3
        )

        genome_end = (
            genome_nt_start
            + (end_aa * 3)
            - 1
        )

        results.append({
            "qseqid": query_id,
            "sseqid": subject_id,
            "label": reference["label"],
            "start_aa": start_aa,
            "end_aa": end_aa,
            "genome_start": genome_start,
            "genome_end": genome_end,
            "pident": row["pident"],
            "bitscore": row["bitscore"],
        })

    return pd.DataFrame(results)


# =============================================================================
# TBLASTN GENOTYPING
# =============================================================================

def run_tblastn(
    query_fasta,
    genotype_reference,
    output_tsv,
    evalue=10,
    min_query_coverage=100,
):
    """
    TBLASTN peptide queries against nucleotide genotype references.

    Up to 50 targets are retained so that multiple near-equivalent genotype
    references can be evaluated using the 98% bitscore rule.
    """

    command = [
        "tblastn",

        "-query",
        str(query_fasta),

        "-subject",
        str(genotype_reference),

        "-out",
        str(output_tsv),

        "-evalue",
        str(evalue),

        "-max_target_seqs",
        "50",

        "-max_hsps",
        "1",

        "-qcov_hsp_perc",
        str(min_query_coverage),

        "-outfmt",
        (
            "6 "
            "qseqid "
            "sseqid "
            "pident "
            "length "
            "qlen "
            "qcovs "
            "bitscore"
        ),
    ]

    run_command(command)


def determine_genotypes(
    query_fasta,
    genotype_reference,
    evalue=10,
    min_query_coverage=100,
    ambiguity_fraction=0.98,
):
    """
    Assign genotype using the 98% top-bitscore rule.

    For each peptide:

        top_score = maximum bitscore

        qualifying references =
            all references with bitscore >=
            ambiguity_fraction * top_score

    Returns a dictionary keyed by peptide ID.

    Each value contains:

        call
        references
        n_references
        top_bitscore
    """

    with tempfile.TemporaryDirectory() as tmpdir:

        blast_output = Path(tmpdir) / "genotype.tsv"

        run_tblastn(
            query_fasta=query_fasta,
            genotype_reference=genotype_reference,
            output_tsv=blast_output,
            evalue=evalue,
            min_query_coverage=min_query_coverage,
        )

        if (
            not blast_output.exists()
            or blast_output.stat().st_size == 0
        ):

            print(
                "\nNo genotype BLAST hits found."
            )

            return {}

        hits = pd.read_csv(
            blast_output,
            sep="\t",
            header=None,
            names=[
                "qseqid",
                "sseqid",
                "pident",
                "length",
                "qlen",
                "qcovs",
                "bitscore",
            ],
        )

    genotype_map = {}

    for query_id, group in hits.groupby(
        "qseqid"
    ):

        query_id = str(
            query_id
        )

        # Highest scoring genotype/reference
        top_score = group[
            "bitscore"
        ].max()

        # Retain all references within 98% of the best score
        qualifying = group[
            group["bitscore"]
            >= ambiguity_fraction * top_score
        ].copy()

        references = sorted(
            qualifying["sseqid"]
            .astype(str)
            .unique()
        )

        if len(references) == 0:

            call = "NoGenotype"

        elif len(references) == 1:

            call = references[0]

        else:

            # Pipe-separated, matching Emma's script. Downstream parsers take
            # the strain as everything between the segment and the coordinate
            # field, so an embedded pipe is expected. A comma would also end up
            # inside a quoted field in the CSV output, which breaks naive readers.
            call = (
                f"Ambiguous{len(references)}:"
                + "|".join(references)
            )

        genotype_map[query_id] = {
            "call": call,
            "references": references,
            "n_references": len(references),
            "top_bitscore": top_score,
        }

    return genotype_map


# =============================================================================
# CPM PEPTIDE COLUMN
# =============================================================================

def find_peptide_column(cpm, requested=None):
    """
    Identify the peptide identifier column.

    Prefer exact `peptide`.

    Also accepts common alternatives used before the final output is made.
    """

    if requested:
        if requested not in cpm.columns:
            raise ValueError(
                f"--cpm-name-column {requested!r} is not in the matrix. "
                f"Available: {list(cpm.columns[:8])}..."
            )
        return requested

    # Matched case-insensitively, most name-like first. `peptide_id` is LAST and
    # deliberately so: in a phip-flow matrix it is a bare row number, not a name,
    # and picking it silently relabels nothing.
    lookup = {str(c).strip().lower(): c for c in cpm.columns}

    for candidate in [
        "peptide",
        "full_name",
        "fullname",
        "peptide_name",
        "name",
        "id",
        "peptide_id",
    ]:

        if candidate in lookup:
            column = lookup[candidate]

            values = cpm[column].astype(str).head(200)
            if values.str.fullmatch(r"\d+").all():
                raise ValueError(
                    f"Column {column!r} holds bare row numbers, not peptide names, "
                    "so nothing would match the library. Phip-flow matrices key "
                    "rows this way and carry the real name in `full_name`. Pass "
                    "--cpm-name-column full_name --cpm-drop-columns peptide_id, or "
                    "prepare the matrix with annotate_plates.py first."
                )

            return column

    raise ValueError(
        "Could not find the peptide identifier column in the CPM matrix. "
        f"Expected `peptide` or a recognised alternative; found {list(cpm.columns[:8])}"
    )


# =============================================================================
# CPM RELABELLING
# =============================================================================

def relabel_cpm(
    cpm_file,
    output_file,
    coordinate_df,
    genotype_map,
    virus_name,
    selected_ids=None,
    name_column=None,
    drop_columns=(),
    segment=None,
    drop_unselected=False,
    coord_map_supplied=False,
):
    """
    Relabel CPM peptide identifiers.

    Final format:

        original|virus|segment|genotype|start-end

    Example:

        pep_00123|cchfv|m|DQ211632_CCHF_IVa_Asia1_M|12345-12510

    Peptides that were selected but could not be placed by BLAST keep the same
    shape, with NOCOORD in the coordinate field:

        pep_00124|cchfv|m|NoGenotype|NOCOORD

    Peptides that were never selected are left with their ORIGINAL name, so the
    matrix stays usable for other viruses. Pass drop_unselected=True to write a
    virus-only matrix instead.
    """

    # The CPM matrix may be comma- or tab-separated. Pick the delimiter from the
    # file extension, as .csv matrices are what the rest of the pipeline produces.
    in_delim = (
        "\t"
        if str(cpm_file).lower().endswith((".tsv", ".txt"))
        else ","
    )

    cpm = pd.read_csv(
        cpm_file,
        sep=in_delim
    )

    peptide_column = find_peptide_column(
        cpm,
        name_column
    )

    # A leftover bare row-index column (phip-flow's `peptide_id`) would be read as
    # a SAMPLE by the rank/call pipeline, with values running 0..n-1 — enormous
    # fake signal that corrupts the mock statistics. Drop it without being asked.
    auto_drop = []
    for column in cpm.columns:
        if column == peptide_column:
            continue
        values = cpm[column].astype(str)
        if values.str.fullmatch(r"\d+").all() and \
                values.astype("int64").tolist() == list(range(len(cpm))):
            auto_drop.append(column)

    drop = [c for c in list(drop_columns) + auto_drop
            if c and c != peptide_column]
    drop = list(dict.fromkeys(drop))
    missing = [c for c in drop if c not in cpm.columns]
    if missing:
        raise ValueError(f"--cpm-drop-columns not in the matrix: {missing}")
    if drop:
        cpm = cpm.drop(columns=drop)
        print(f"  dropped columns: {', '.join(drop)}")

    # -------------------------------------------------------------------------
    # Coordinate lookup
    # -------------------------------------------------------------------------

    coordinate_lookup = {}

    for _, row in coordinate_df.iterrows():

        coordinate_lookup[
            str(row["qseqid"])
        ] = row

    # -------------------------------------------------------------------------
    # Relabel each CPM peptide
    # -------------------------------------------------------------------------

    new_names = []
    keep_row = []

    n_relabelled = 0
    n_placed = 0
    n_unselected = 0

    for peptide in cpm[
        peptide_column
    ].astype(str):

        # ---------------------------------------------------------------------
        # Peptides this run never selected keep their original name.
        # ---------------------------------------------------------------------
        if (
            selected_ids is not None
            and peptide not in selected_ids
        ):

            new_names.append(peptide)
            keep_row.append(not drop_unselected)
            n_unselected += 1

            continue

        coordinate = coordinate_lookup.get(
            peptide
        )

        # ---------------------------------------------------------------------
        # Genotype
        # ---------------------------------------------------------------------

        genotype_info = genotype_map.get(
            peptide
        )

        genotype = (
            genotype_info["call"]
            if genotype_info is not None
            else "NoGenotype"
        )

        # ---------------------------------------------------------------------
        # Segment / protein field
        # ---------------------------------------------------------------------

        if segment:

            segment_field = segment

        elif coordinate is not None:

            segment_field = str(
                coordinate["label"]
            )

        else:

            segment_field = None

        # ---------------------------------------------------------------------
        # Coordinates
        # ---------------------------------------------------------------------

        coordinate_field = "NOCOORD"

        if (
            coordinate is not None
            and segment_field != "NOCOORD"
        ):

            genome_start = coordinate["genome_start"]
            genome_end = coordinate["genome_end"]

            # Without --coord-map there are no genome nucleotide coordinates, but
            # the BLAST amino-acid positions on the coordinate reference are still
            # known. Fall back to those rather than discarding the placement.
            #
            # Only when NO coord map was given, though. If one was supplied and
            # this subject is missing from it, falling back would silently mix
            # amino-acid and nucleotide coordinates in the same column.
            if (
                not coord_map_supplied
                and (
                    pd.isna(genome_start)
                    or pd.isna(genome_end)
                )
            ):

                genome_start = coordinate["start_aa"]
                genome_end = coordinate["end_aa"]

            if not (
                pd.isna(genome_start)
                or pd.isna(genome_end)
            ):

                coordinate_field = (
                    f"{int(genome_start)}-"
                    f"{int(genome_end)}"
                )

                n_placed += 1

        # ---------------------------------------------------------------------
        # Assemble
        # ---------------------------------------------------------------------

        if segment_field is None or segment_field == "NOCOORD":

            # Nothing to put in the segment slot: keep the short legacy shape.
            new_names.append(
                f"{peptide}|{virus_name}|NOCOORD"
            )

        else:

            new_names.append(
                f"{peptide}|"
                f"{virus_name}|"
                f"{segment_field}|"
                f"{genotype}|"
                f"{coordinate_field}"
            )

        keep_row.append(True)
        n_relabelled += 1

    if n_relabelled == 0:
        example_matrix = cpm[peptide_column].astype(str).head(2).tolist()
        example_ids = list(selected_ids)[:2] if selected_ids else []
        raise ValueError(
            "No peptide in the CPM matrix matched the library.\n"
            f"  matrix {peptide_column!r} looks like: {example_matrix}\n"
            f"  library IDs look like:              {example_ids}\n"
            "The two must be the same identifier. Check --cpm-name-column."
        )

    print(
        f"\n  rows read        {len(new_names):,}"
        f"\n  relabelled       {n_relabelled:,}"
        f"\n    of which placed {n_placed:,}"
        f"\n  left unchanged   {n_unselected:,}"
        + (
            f"  (dropped from output)"
            if drop_unselected
            else ""
        )
    )

    # -------------------------------------------------------------------------
    # Final column name expected by downstream pipeline
    # -------------------------------------------------------------------------

    cpm[
        peptide_column
    ] = new_names

    if peptide_column != "peptide":

        cpm = cpm.rename(
            columns={
                peptide_column: "peptide"
            }
        )

    # The peptide column must come FIRST. Downstream readers use index_col=0, and
    # the name column is not always first in the input (phip-flow matrices carry
    # it last, in `full_name`).
    if list(cpm.columns)[0] != "peptide":

        cpm = cpm[
            ["peptide"]
            + [c for c in cpm.columns if c != "peptide"]
        ]

    # -------------------------------------------------------------------------
    # Drop unselected rows if a virus-only matrix was requested
    # -------------------------------------------------------------------------

    if drop_unselected:

        cpm = cpm[
            pd.Series(keep_row, index=cpm.index)
        ]

    # -------------------------------------------------------------------------
    # Write output
    # -------------------------------------------------------------------------

    out_delim = (
        "\t"
        if str(output_file).lower().endswith((".tsv", ".txt"))
        else ","
    )

    cpm.to_csv(
        output_file,
        sep=out_delim,
        index=False
    )

    print(
        f"\nWrote relabelled CPM:"
        f"\n  {output_file}"
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Select virus peptides from a translated PhIP-seq library, "
            "map them to a coordinate reference, optionally genotype them, "
            "and relabel a CPM matrix."
        )
    )

    # =========================================================================
    # INPUTS
    # =========================================================================

    parser.add_argument(
        "--library",
        required=True,
        help=(
            "Translated peptide FASTA library. "
            "Use the actual translated peptide library used for PhIP-seq."
        ),
    )

    parser.add_argument(
        "--coord-ref",
        required=True,
        help=(
            "Protein FASTA used as the coordinate reference."
        ),
    )

    parser.add_argument(
        "--coord-map",
        required=False,
        help=(
            "Tab-separated coordinate map containing: "
            "subject_id, label, genome_nt_start."
        ),
    )

    parser.add_argument(
        "--geno-ref",
        required=False,
        help=(
            "Optional nucleotide FASTA containing genotype/reference "
            "sequences for TBLASTN."
        ),
    )

    parser.add_argument(
        "--cpm",
        required=True,
        help="CPM matrix to relabel.",
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output relabelled CPM matrix.",
    )

    parser.add_argument(
        "--cpm-name-column",
        required=False,
        default=None,
        help=(
            "Column of the CPM matrix holding the peptide NAME, when it is not "
            "the first column or is not auto-detected. Phip-flow matrices often "
            "key rows by a bare row number in `peptide_id` and carry the real "
            "name in `full_name`."
        ),
    )

    parser.add_argument(
        "--cpm-drop-columns",
        required=False,
        default="",
        help=(
            "Comma-separated columns to remove before writing. Anything left in "
            "that is not the peptide column is treated as a SAMPLE by the "
            "rank/call pipeline, so a row-index column such as `peptide_id` must "
            "be dropped or it becomes a fake sample with enormous signal."
        ),
    )

    parser.add_argument(
        "--virus",
        required=True,
        help="Virus name used in the final peptide identifiers.",
    )

    parser.add_argument(
        "--segment",
        required=False,
        default=None,
        help=(
            "Optional segment or protein label placed in the third field of the "
            "peptide name, e.g. m or s. Downstream plotting expects "
            "virus|segment|genotype|coords. If omitted, the BLAST subject ID is "
            "used for placed peptides."
        ),
    )

    parser.add_argument(
        "--drop-unselected",
        action="store_true",
        help=(
            "Write a virus-only matrix, dropping peptides this run did not "
            "select. By default those rows are kept with their original names "
            "so the matrix stays usable for other viruses."
        ),
    )

    # =========================================================================
    # PEPTIDE SELECTION
    # =========================================================================

    parser.add_argument(
        "--search-term",
        required=False,
        help=(
            "Regular expression used to select peptides by FASTA header. "
            "This is the DEFAULT selection method. "
            "Use | for alternatives."
        ),
    )

    parser.add_argument(
        "--sequence-search",
        action="store_true",
        help=(
            "BLAST the entire translated library against --coord-ref "
            "to identify candidate virus peptides. "
            "Use for poorly named/inconsistently annotated libraries."
        ),
    )

    parser.add_argument(
        "--min-identity",
        type=float,
        default=None,
        help=(
            "Optional minimum percent identity for --sequence-search."
        ),
    )

    # =========================================================================
    # BLAST PARAMETERS
    # =========================================================================

    parser.add_argument(
        "--evalue",
        type=float,
        default=10,
        help="BLAST e-value cutoff. Default: 10.",
    )

    parser.add_argument(
        "--min-query-coverage",
        type=float,
        default=100,
        help=(
            "Minimum BLAST query coverage. "
            "Default: 100."
        ),
    )

    parser.add_argument(
        "--ambiguity-fraction",
        type=float,
        default=0.98,
        help=(
            "Retain genotype references within this fraction of the "
            "top bitscore. Default: 0.98."
        ),
    )

    args = parser.parse_args()

    # =========================================================================
    # PATHS
    # =========================================================================

    library = Path(
        args.library
    )

    coord_ref = Path(
        args.coord_ref
    )

    cpm = Path(
        args.cpm
    )

    output = Path(
        args.output
    )

    # =========================================================================
    # CHECK SELECTION MODE
    # =========================================================================

    if (
        args.sequence_search
        and args.search_term
    ):

        raise ValueError(
            "Use either --search-term OR --sequence-search, "
            "not both."
        )

    if (
        not args.sequence_search
        and not args.search_term
    ):

        raise ValueError(
            "You must provide --search-term unless "
            "--sequence-search is being used."
        )

    # =========================================================================
    # TEMPORARY WORKSPACE
    # =========================================================================

    with tempfile.TemporaryDirectory() as tmpdir:

        tmpdir = Path(
            tmpdir
        )

        selected_fasta = (
            tmpdir / "virus_peptides.fa"
        )

        coordinate_blast = (
            tmpdir / "coordinate_blast.tsv"
        )

        # =====================================================================
        # STEP 1: SELECT PEPTIDES
        # =====================================================================

        if args.sequence_search:

            print(
                "\n"
                "============================================================\n"
                "PEPTIDE SELECTION: SEQUENCE SEARCH\n"
                "============================================================"
            )

            print(
                "The entire translated library will be BLASTed "
                "against the target reference."
            )

            selected = sequence_search_library(
                library_fasta=library,
                reference_fasta=coord_ref,
                output_fasta=selected_fasta,
                evalue=args.evalue,
                min_query_coverage=args.min_query_coverage,
                min_identity=args.min_identity,
            )

        else:

            print(
                "\n"
                "============================================================\n"
                "PEPTIDE SELECTION: NAME / TERM SEARCH\n"
                "============================================================"
            )

            print(
                f"Search term: {args.search_term}"
            )

            library_records = read_fasta(
                library
            )

            selected = select_by_name(
                library_records,
                args.search_term
            )

            write_fasta(
                selected,
                selected_fasta
            )

            print(
                f"Selected "
                f"{len(selected):,} / "
                f"{len(library_records):,} "
                f"library peptides."
            )

        if not selected:

            raise RuntimeError(
                "No target-virus peptides were selected.\n\n"
                "If the library is correctly named, check --search-term.\n"
                "If the library is poorly named, try --sequence-search."
            )

        # =====================================================================
        # STEP 2: COORDINATE BLAST
        # =====================================================================

        print(
            "\n"
            "============================================================\n"
            "COORDINATE BLASTP\n"
            "============================================================"
        )

        run_blastp(
            query_fasta=selected_fasta,
            reference_fasta=coord_ref,
            output_tsv=coordinate_blast,
            evalue=args.evalue,
            min_query_coverage=args.min_query_coverage,
            max_target_seqs=1,
        )

        if (
            not coordinate_blast.exists()
            or coordinate_blast.stat().st_size == 0
        ):

            raise RuntimeError(
                "No coordinate BLAST hits were found."
            )

        # Keep the coordinate BLAST next to the output. It is the only record of
        # the per-peptide identity distribution, and choosing --min-identity for
        # a new virus means reading that distribution: PROTOCOL_relabelling_any
        # _library.md documents exactly that procedure, against exactly this
        # file. Until now it lived in a TemporaryDirectory and was deleted on
        # exit, so the one step the docs call out as the thing people get wrong
        # could not be performed at all.
        coords_kept = Path(args.output).with_suffix("").as_posix() + "_peptide_coords.tsv"
        # Written WITH a header, unlike BLAST's own outfmt 6, because the
        # procedure in PROTOCOL step 2 reads it with csv.DictReader and asks for
        # the 'pident' column by name.
        coords_header = ("qseqid\tsseqid\tsstart\tsend\tpident\t"
                         "length\tqlen\tqcovs\tbitscore\n")
        try:
            with open(coords_kept, "w") as kept:
                kept.write(coords_header)
                with open(coordinate_blast) as raw:
                    shutil.copyfileobj(raw, kept)
            print(f"\nCoordinate BLAST kept at: {coords_kept}")
            print("  Use it to choose --min-identity for a new virus; see "
                  "PROTOCOL_relabelling_any_library.md, step 2.")
        except OSError as exc:
            print(f"\nWarning: could not keep the coordinate BLAST ({exc}).")

        coordinate_hits = pd.read_csv(
            coordinate_blast,
            sep="\t",
            header=None,
            names=[
                "qseqid",
                "sseqid",
                "sstart",
                "send",
                "pident",
                "length",
                "qlen",
                "qcovs",
                "bitscore",
            ],
        )

        print(
            f"Coordinate BLAST produced "
            f"{len(coordinate_hits):,} hits."
        )

        # =====================================================================
        # STEP 3: COORDINATE MAPPING
        # =====================================================================

        print(
            "\n"
            "============================================================\n"
            "COORDINATE MAPPING\n"
            "============================================================"
        )

        if args.coord_map:

            coord_map = load_coord_map(
                args.coord_map
            )

            coordinate_df = map_coordinates(
                coordinate_hits,
                coord_map
            )

        else:

            print(
                "No --coord-map supplied."
            )

            coordinate_df = (
                coordinate_hits.copy()
            )

            coordinate_df["label"] = (
                coordinate_df["sseqid"]
            )

            coordinate_df["start_aa"] = (
                coordinate_df[
                    ["sstart", "send"]
                ].min(axis=1)
            )

            coordinate_df["end_aa"] = (
                coordinate_df[
                    ["sstart", "send"]
                ].max(axis=1)
            )

            coordinate_df[
                "genome_start"
            ] = pd.NA

            coordinate_df[
                "genome_end"
            ] = pd.NA

        # =====================================================================
        # STEP 4: GENOTYPING
        # =====================================================================

        if args.geno_ref:

            print(
                "\n"
                "============================================================\n"
                "GENOTYPING\n"
                "============================================================"
            )

            print(
                "Using 98% of top bitscore ambiguity rule."
            )

            genotype_map = determine_genotypes(
                query_fasta=selected_fasta,
                genotype_reference=args.geno_ref,
                evalue=args.evalue,
                min_query_coverage=args.min_query_coverage,
                ambiguity_fraction=args.ambiguity_fraction,
            )

            print(
                f"Genotype calls generated for "
                f"{len(genotype_map):,} peptides."
            )

        else:

            print(
                "\nNo --geno-ref supplied."
            )

            genotype_map = {}

        # =====================================================================
        # STEP 5: RELABEL CPM
        # =====================================================================

        print(
            "\n"
            "============================================================\n"
            "RELABEL CPM\n"
            "============================================================"
        )

        relabel_cpm(
            cpm_file=cpm,
            output_file=output,
            coordinate_df=coordinate_df,
            genotype_map=genotype_map,
            virus_name=args.virus,
            selected_ids={
                record["id"]
                for record in selected
            },
            name_column=args.cpm_name_column,
            drop_columns=[c.strip() for c in args.cpm_drop_columns.split(",")
                          if c.strip()],
            segment=args.segment,
            drop_unselected=args.drop_unselected,
            coord_map_supplied=bool(args.coord_map),
        )

    print(
        "\n============================================================"
    )

    print(
        "DONE"
    )

    print(
        "============================================================"
    )


if __name__ == "__main__":
    main()