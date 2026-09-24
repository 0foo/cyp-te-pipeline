# RECONSTRUCTED — not a lab original.
#
# NEW_Step_5_Replace_gff_Names_with_Dmelanogaster_1_9.py and
# Step_2_Remove_Duplicate_gene_names_10_31.py both `from config import ...`, but the lab's
# config.py was never kept. This file supplies the names they import, written on 2026-09-24.
#
# Validated: with the _10_31 HOG table below, NEW_Step_5 run on the Zenodo D. ananassae
# annotation and passed through repeatOpp.py reproduces the lab's committed filtered.gff
# exactly (166 Cyp gene lines). See docs/01-data-preparation.md.
#
# Every path can be overridden with an environment variable of the same name.

import os

HERE = os.path.dirname(os.path.abspath(__file__))

# The HOG table: one row per orthogroup, one column per species, genes comma-separated.
# Use the _10_31 version. The _1_9 version has no line breaks (every row run together) and
# cannot be parsed; the undated version does not reproduce the lab's results.
STEP5_INPUT = os.environ.get(
    "STEP5_INPUT",
    os.path.join(
        HERE, "..", "..", "datasets", "hog-table",
        "HOG_OG_association_gene_names_without_duplicates_10_31.tsv",
    ),
)

# Folder holding the Zenodo annotations, one <SPECIES>_final.gff per species (unzipped).
# NEW_Step_5 processes every species whose file is present here, so put only the species
# you want in this folder.
GFF_INPUT_DIR = os.environ.get("GFF_INPUT_DIR", os.path.join(HERE, "input"))

# Where <SPECIES>_final_withDmelNames.gff files are written.
GFF_OUTPUT_DIR = os.environ.get("GFF_OUTPUT_DIR", os.path.join(HERE, "output"))

# Step 2 only (not part of the happy path; see docs/01-data-preparation.md).
STEP2_INPUT = os.environ.get("STEP2_INPUT", "HOG_OG_association_gene_names.tsv")
STEP2_OUTPUT = os.environ.get("STEP2_OUTPUT", "HOG_OG_association_gene_names_without_duplicates.tsv")


def _dmel_column_index(path):
    with open(path) as f:
        return f.readline().rstrip("\n").split("\t").index("DROSOPHILA_MELANOGASTER")


# Column of the HOG table holding the D. melanogaster gene symbols.
DMEL_COLUMN_INDEX = _dmel_column_index(STEP5_INPUT)


def get_gff_input_path(label, step=None):
    return os.path.join(GFF_INPUT_DIR, f"{label}_final.gff")


def get_gff_output_path(label):
    os.makedirs(GFF_OUTPUT_DIR, exist_ok=True)
    return os.path.join(GFF_OUTPUT_DIR, f"{label}_final_withDmelNames.gff")


def clean_field(x):
    # Multi-gene cells in the HOG table are written as "gene-A, gene-B" (quoted, with a
    # space after each comma). Stripping whitespace and quotes recovers each gene name.
    # Stripping whitespace is also what the lab's version must have done: it explains why
    # the _1_9 table, written by Step 2 through clean_field, lost every line break.
    return x.strip().strip('"').strip()
