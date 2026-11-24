# validations

**New (not lab files, 2026-11-23/24).** These scripts check what the pipeline's own scripts acquired:
which Cyp genes they found and which TEs they counted. They read the outputs of a full pipeline run.
They don't change or re-run any lab script, and nothing in `data-preparation/` or `data-analysis/` depends on them.
They produced the Google Sheet "Cyp-TE validation 2026-11-24".

Paths are hard-coded to the run in `/bioinformatics/atallah-repeats` on mother:

- `rmodeler/out/` (stage 1)
- `3-pair-te-with-cyp/` (stages 2–3)
- `4-combine-and-compare/` (stage 4)
- `summary/`, where the outputs are written

Run with `python3 -B`, so that no `__pycache__` is written into the repo.

| Script | What it does |
|---|---|
| `make-summary.py` | Per-species and per-Cyp-locus tables from stages 1–3. Writes `summary/per-species.tsv`, `cyp-genes.tsv` and `xenobiotic-genes.tsv`; `make-validation-sheet.py` reads the first two. Standard library only. |
| `make-validation-sheet.py` | Builds `summary/cyp-te-validation.xlsx`, one row per species (301). **Cyp tab:** Cyp loci `repeatOpp.py` selected; unique names after `CleanAnnasse.py`; non-P450 loci (cyclophilins, Cpr...); % absent against 84 P450s per genome (Good et al. 2014: 74–94 across 12 *Drosophila*); the 4 stable Halloween genes and 7 conserved single-copy Cyps with known conserved functions, both as found by the scripts and as present in the renamed annotation. The gene sets come from Good et al. 2014, *GBE* 6:1118 ([doi:10.1093/gbe/evu083](https://doi.org/10.1093/gbe/evu083)): genes that are single-copy in all 12 *Drosophila* genomes (Table 1) and to which the paper assigns a function, plus *shadow*/Cyp315a1 (Warren et al. 2002). See the comment above `CONSERVED` in the script. **Transposons tab:** genome-wide TE insertions and coverage from RepeatMasker `.out`; % of hits that are not true TEs; repeat-library size; true TEs near Cyps (stage 3); the "TEs" stage 4 counted and the share of them that are not true TEs; TE density near Cyps divided by genome-wide density. Needs `openpyxl`, e.g. `nix-shell -p 'python3.withPackages(p: [p.openpyxl])'`. Reads about 6 GB of `.out` files, roughly 10 minutes. |
| `column_labels.py` | Descriptive column headers (`LABELS`) and the columns shown on the Transposons tab (`TE_KEEP`). Edit this file to rename or show columns. |
| `relabel.py` | Re-applies `column_labels.py` to an already-built `cyp-te-validation.xlsx`, without re-reading the RepeatMasker files. |

Each data tab in the workbook is followed by a "– columns" tab that defines every column.

## Where the transposon numbers come from

Every TE number traces back to one RepeatMasker run per species. Three different counts are then taken from it:

1. **Genome-wide** (Transposons tab, "Genome-wide …" and "% of genome covered by TEs").
   - Stage 1: `data-preparation/1-repeat-library-and-masking/worker.sh`, already in the repo, runs
     `BuildDatabase`, then `RepeatModeler -database <sample>`, then `RepeatMasker -lib <sample>-families.fa -xsmall` in the Dfam TE Tools container.
     The exact commands are in `rmodeler/logs/*.log` (RepeatModeler 2.0.9, RepeatMasker 4.2.4).
   - Output: `<sample>.rm.out`, `.rm.tbl` and `-families.fa`, all from the Drive archive `rmodeler.tar.gz`.
   - `pipeline-run/1-repeat-library-and-masking/check-rm-out.py` and `grade-rm-out.py` verified all 319 `.rm.out` files before use.
   - `make-validation-sheet.py` then reads each `.rm.out` row and classifies it by column 11 (class/family).
     Simple_repeat, Low_complexity, Satellite, rRNA/tRNA/snRNA/scRNA/srpRNA and ARTEFACT are counted as **not true TEs**; everything else is a TE.
     - "TE insertions" counts distinct element IDs (column 15), so fragments of one insertion count once.
     - TE bp is the sum of the TE row lengths.
2. **Near Cyps, stage 3** ("Near Cyps (stage 3) …").
   - `pipeline-run/3-pair-te-with-cyp/make-pairs.py` pairs each RepeatMasker sample with its annotation and checks they come from the same assembly.
   - `run-pairing.sh` then runs the lab's `repeatOpp.py` → `Locate_TE.py` → `CleanAnnasse.py`, unchanged, as path-patched copies.
     Together they keep every RepeatMasker row lying fully inside a Cyp gene ±3 kb.
   - `make-summary.py` removes repeatOpp's duplicate rows (one count per locus × RepeatMasker row) and splits the hits into TE and not-TE using the same class rule.
3. **Stage 4** ("Stage 4 …").
   - `pipeline-run/4-combine-and-compare/run-build.sh` runs `build_tfbs_te_gff.py` unchanged. FIMO runs in the `fimo-image/Dockerfile` image.
     Every stage 3 hit becomes a `mobile_genetic_element` row tagged "Within range of <gene>".
   - `run-compare.sh` runs the three `compare_te_cyp_*.py` scripts unchanged. They count every such row as a "TE", simple repeats included.
   - `make-validation-sheet.py` re-counts those rows with the same rule as `compare_te_cyp_exposure.py`. The totals match the scripts' output for all 235 tested species.
     It then splits them by `repeat_class`.
4. **Density ratio.** (TE hits near Cyps per 10 kb of Cyp window) ÷ (TE hits per 10 kb of genome), both counted in RepeatMasker rows.

## Other folders

- `pipeline-run/`: the runners and checks used for the full run, copied from `/bioinformatics/atallah-repeats/<stage>/`.
  - `2-rename-genes/run-rename.sh`: runs NEW_Step_5 in batches.
  - `4-combine-and-compare/make-exposure.py`: high/low exposure groups from the ecology sheet.
  - `4-combine-and-compare/check-genome-vs-gff.sh`: genome-to-annotation check.

  These only pass paths and flags to the lab scripts, which they never edit.

  On 2026-11-24 the header comments of these repo copies were expanded (comments only), so they no
  longer match the copies in `/bioinformatics/atallah-repeats/` byte for byte. The exception is
  `run-pairing.sh`: its sed line numbers were updated for the lab scripts' new documentation headers
  (46/48, 49/52, 43; they were 4/6, 2/5, 3), and the repo copy was copied back over the live one so the
  two stay identical.
- `sheets/`: `make-sheet.py` builds "Cyp-TE pipeline results" (per-species data plus stage 4 statistics), and
  `make-hits-sheet.py` builds "Cyp-TE per-species hits and stage 4". Both need `openpyxl`.

Findings at the time of writing (see the sheet for per-species values):

- The scripts find a median of 15 unique Cyp names per species, about 82% fewer than expected.
- None of the loci the scripts found are non-P450s.
- The scripts find no Halloween genes in any species except *D. melanogaster*. The lab gene list names them by number (Cyp306a1...), but the renamed annotations name them by symbol (phm...).
- About two-thirds of the hits stage 4 counts as "TEs" are microsatellite or low-complexity repeats.
- Cyp regions hold about 0.57× the genome-wide TE density.
