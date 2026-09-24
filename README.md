# cyp-te-pipeline

The happy path through the Atallah lab's Cyp / transposable-element pipeline: which
transposable elements (TEs) sit in or near cytochrome P450 (Cyp) genes, species by species,
and whether that differs between species with high and low xenobiotic exposure.

Every script here is a **byte-for-byte copy** of a lab original (checksums in
[`docs/03-provenance.md`](docs/03-provenance.md)). The one exception is
`data-preparation/3-rename-genes/config.py`, a reconstruction of a lost lab file, labelled as
such and validated against the lab's own results. The originals stay untouched in
`~/projects/cyp-project-forensics` and `~/projects/repeat-modeler-automation`.

## The pipeline

```
                    DATA PREPARATION                                 DATA ANALYSIS

genome FASTA ──► 1-repeat-library-and-masking ──► <species>.rm.out ─┐
                 (RepeatModeler + RepeatMasker)                     │
                                                                    ├─► 1-pair-te-with-cyp ──► D_<species>GenesAffectedByTE.txt
Zenodo <SPECIES>_final.gff ─┐                                       │   (repeatOpp → Locate_TE          │
                            ├─► 3-rename-genes ──► <SPECIES>_final_ ┘    → CleanAnnasse)               │
HOG table (_10_31) ─────────┘   (NEW_Step_5)      withDmelNames.gff ─────────────┐                      │
     ▲                                                                           ▼                      ▼
     └── 2-orthogroup-table (Step_1, Step_2) — reference only           2-combine-and-compare ◄─────────┘
                                                                         build_tfbs_te_gff → compare_te_cyp_*
                                                                                   │
                                                                         CSV + report + plot
```

## Layout

| Folder | Contents |
|---|---|
| [`data-preparation/1-repeat-library-and-masking/`](data-preparation/1-repeat-library-and-masking/) | `worker.sh`, `rm-manager.sh`, `rmodeler.conf.example` — genome → repeat library → TE locations (`.rm.out`) |
| [`data-preparation/2-orthogroup-table/`](data-preparation/2-orthogroup-table/) | `Step_1_…`, `Step_2_…` — how the HOG table was built. **Reference only**: they cannot be re-run, and the table they made is the input to step 3 |
| [`data-preparation/3-rename-genes/`](data-preparation/3-rename-genes/) | `NEW_Step_5_…` + reconstructed `config.py` — relabels each species' genes with *D. melanogaster* ortholog names |
| [`data-analysis/1-pair-te-with-cyp/`](data-analysis/1-pair-te-with-cyp/) | `repeatOpp.py`, `Locate_TE.py`, `CleanAnnasse.py`, `Reg_Gene_Full.txt` — Cyp genes × nearby TEs |
| [`data-analysis/2-combine-and-compare/`](data-analysis/2-combine-and-compare/) | `build_tfbs_te_gff.py`, `compare_te_cyp_{exposure,cncc,xenobiotic}.py` — combined annotation, then the cross-species statistics |
| [`docs/`](docs/) | how to run each step, and where every file came from |

## Start here

1. [`docs/01-data-preparation.md`](docs/01-data-preparation.md) — getting a species from raw
   downloads to a TE table and a renamed annotation.
2. [`docs/02-data-analysis.md`](docs/02-data-analysis.md) — from those two files to the
   cross-species comparison.
3. [`docs/03-provenance.md`](docs/03-provenance.md) — what was copied from where, what was
   reconstructed, how it was validated, and why this renaming script and not the other one.

## Why this is the "correct" path

The lab had two gene-renaming scripts. This pipeline uses Ayush's
`NEW_Step_5_Replace_gff_Names_with_Dmelanogaster_1_9.py`, not Duy's `ReVamp_Final.py`,
because ReVamp mis-parses orthogroup cells that list several genes and silently leaves those
genes unrenamed. On *D. ananassae* that drops 34 of 91 Cyp loci; across the lab's finished
tables, 19 species were built with ReVamp and are missing such loci. Details:
[`docs/03-provenance.md`](docs/03-provenance.md).

## Requirements

- Python 3. **pandas 2.x** for `repeatOpp.py` (imported, never used); **matplotlib** for the
  comparison plots (skippable with `--no-plot`). Everything else is standard library.
- Docker and the `dfam/tetools` image (step 1 of data preparation)
- MEME Suite's `fimo`, locally or via Docker (for transcription-factor scanning in
  `build_tfbs_te_gff.py`; skippable with `--skip-tfbs`)
- Local disk for all work — not network mounts
