# Notes — open gaps

What the repository still lacks, as of 2026-09-24. Everything needed for the validated part
(renaming *D. ananassae* and pairing its Cyp genes with TEs, checked against the lab's results)
is already here. Items marked *new file* would be written fresh, not copied from the lab, and
labelled as such, like `data-preparation/3-rename-genes/config.py`.

## Would block a run

1. ~~**Where genomes come from.**~~ *Recorded 2026-09-24:* `docs/01-data-preparation.md` now
   gives the source — Zenodo record 18453526, `genomes.tar.gz` (19.6 GB, md5
   `bca079304da4dbe8e0c9998fc049eb03`), with NCBI GCF_017639315.1 as the single-species
   alternative for *D. ananassae* and `--sequence-source ncbi` to skip the FASTA altogether.
   No longer blocking.
2. **Species config for the comparison** (*new file*). `compare_te_cyp_exposure.py`,
   `compare_te_cyp_cncc.py` and `compare_te_cyp_xenobiotic.py` need an INI file giving each
   species' combined GFF3 and whether it is `high` or `low` exposure. The lab's example
   (`te_cyp_species_config.example.ini`) was never kept, and the exposure assignments are not
   recorded anywhere — only *D. santomea* is noted as low. Needs a template to fill in.
3. **Gene list for the xenobiotic comparison** (*new file*). `--gene-list` takes one gene name
   per line, `#` for comments. The lab's example
   (`xenobiotic_resistance_cyp_genes.example.txt`) was never kept. Needs a template.

## Would help

4. **`requirements.txt`** (*new file*): `pandas<3` (for `repeatOpp.py`; pandas 3 changes string
   handling) and `matplotlib` (comparison plots). Currently only described in the README.
5. **A check script** (*new file*) that runs the *D. ananassae* validation in one command:
   rename with `NEW_Step_5`, then `repeatOpp.py` → `Locate_TE.py` → `CleanAnnasse.py`, and
   compare each output with `datasets/ananassae-reference/`. Today the procedure is only in
   the docs. Not yet run end to end from inside this repository.
6. **Broken links in the copied automation docs.** `docs/repeat-modeler-automation/README.md`
   links to `docs/SOURCE.md` and `SOURCE.md` links to `../README.md` — their original folder
   layout. Both sit in one folder here, so neither link resolves.

## Limitation, not a gap

- **Hardcoded paths in the pairing scripts.** `repeatOpp.py`, `Locate_TE.py` and
  `CleanAnnasse.py` keep the lab's Windows paths and must be edited before each run (lines
  listed in `docs/02-data-analysis.md`). A hands-off option that keeps the lab copies unchanged:
  a small wrapper that writes path-patched copies at run time.

## Deliberately not included

- Annotations for the other 300 species — 1.76 GB, downloadable with a checksum
  (`docs/01-data-preparation.md`).
- RepeatMasker `.out` files for other species — the lab kept only a few; step 1 produces them.
- JASPAR motifs — `build_tfbs_te_gff.py` downloads them at run time, or takes
  `--jaspar-motif-file` for offline use.
