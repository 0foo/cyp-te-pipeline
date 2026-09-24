# datasets

Inputs the pipeline needs that cannot be downloaded anywhere else, plus a complete reference
set for one species so every step can be checked against the lab's own results. All files are
unchanged copies; origins and checksums are in
[`../docs/03-provenance.md`](../docs/03-provenance.md).

## `hog-table/`

`HOG_OG_association_gene_names_without_duplicates_10_31.tsv` (72 MB) — one row per
hierarchical orthogroup (HOG), one column per species, each cell the species' genes in that
HOG (comma-separated; cells with several genes are quoted). The `DROSOPHILA_MELANOGASTER`
column gives each HOG's *D. melanogaster* gene name(s).

Read by `data-preparation/3-rename-genes/` (its `config.py` points here by default). Made by
the lab from OrthoFinder-style orthogroups with the scripts in
`data-preparation/2-orthogroup-table/`, which can no longer be re-run.

## `ananassae-reference/`

Everything for *D. ananassae*, the one species the lab left a full worked example for:

| File | What it is | Pipeline role |
|---|---|---|
| `DROSOPHILA_ANANASSAE_final.gff.gz` | the Zenodo annotation (gzipped, as distributed) | input to renaming |
| `Drosophila_ananassae.GCF_017639315.1.rm.fna.out` | RepeatMasker output: every TE in the genome | input to `Locate_TE.py` |
| `filtered.gff` | the lab's `repeatOpp.py` output — 166 Cyp gene lines | expected result |
| `GenesAffectedByTEs.txt` | the lab's `Locate_TE.py` output — 900 rows | expected result |
| `DAnasse_TE_Cyp.txt` | the lab's `CleanAnnasse.py` output — the finished table | expected result |

The three expected results have Windows (CRLF) line endings; strip them (`tr -d '\r'`) before
comparing.

Not included: the *D. ananassae* genome FASTA (only needed for `build_tfbs_te_gff.py`; NCBI
accession GCF_017639315.1, or use `--sequence-source ncbi`), and annotations for the other 300
species (Zenodo record 18453526, see `docs/01-data-preparation.md`).
