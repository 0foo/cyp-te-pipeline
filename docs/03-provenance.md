# Provenance

Where every file in this folder came from, what was reconstructed, and the evidence behind
choosing this path. Assembled 2026-09-24.

Everything needed is inside this repository. The paths below record where each file was
copied *from* — the lab's archive (`cyp-project-forensics`) and the automation repository
(`repeat-modeler-automation`, also at github.com/0foo/repeat-modeler-automation). Neither was
changed, and neither is needed to use this repository. The full forensic write-up is
`docs/deep/06-final-final-gff.md` in `cyp-project-forensics`.

---

## Copied files

Each script below is a copy of a lab (or automation) original that differs from it **only by a
documentation header** added on 2026-11-24 at the user's request. The header explains how the script
works, its inputs and outputs, the theory behind it and its known defects. Nothing else in any copy
was changed: not a code line, not whitespace, not line endings (the CRLF files carry CRLF headers).
`Reg_Gene_Full.txt` has no header and is still a byte-for-byte copy.

The header is one contiguous block between two marker lines:

```
# ===== DOCUMENTATION HEADER (added 2026-11-24; not part of the lab original) =====
# ...
# ===== END DOCUMENTATION HEADER =====
```

It sits at the top of each file, after the `#!` line where there is one. The exception is
`worker.sh`, where it sits after the script's own header comment (after original line 147),
because `usage()` prints lines 3–34 of the file with `sed` and an insertion above line 34 would
change `--help`.

**Recovering an original.** Delete the block, markers included:

```bash
sed '/^# ===== DOCUMENTATION HEADER (added 2026-11-24/,/^# ===== END DOCUMENTATION HEADER =====/d' FILE
```

**Provenance check.** "SHA-256 (original)" is checked on the file with the header block removed;
the first 12 hex digits must match. This was verified for all 13 scripts on 2026-11-24, both
against these values and against the last header-free commit (`git show fb2f8a1:<path>`).
"SHA-256 (as now, with header)" identifies the current file and changes if the header is edited.

| Here | Original | SHA-256 (original; file with header removed) | SHA-256 (as now, with header) |
|---|---|---|---|
| `data-preparation/1-repeat-library-and-masking/worker.sh` | `repeat-modeler-automation/worker.sh` | `bbbe404c8c80` | `d9d9780858a7` |
| `data-preparation/1-repeat-library-and-masking/rm-manager.sh` | `repeat-modeler-automation/rm-manager.sh` | `b61435b5e58b` | `20ff98b07ec7` |
| `data-preparation/1-repeat-library-and-masking/rmodeler.conf.example` | `repeat-modeler-automation/rmodeler.conf.example` | `8d551841c8e0` | `31743c1b0573` |
| `data-preparation/2-orthogroup-table/Step_1_Replace_RNA_identifiers_with_gene_names_Redo.py` | `cyp-project-forensics/to_organize/` (same name) | `74d53b16c85f` | `904bbacb123f` |
| `data-preparation/2-orthogroup-table/Step_2_Remove_Duplicate_gene_names_10_31.py` | `cyp-project-forensics/to_organize/` (same name) | `1d06b9a8d914` | `9062205f6c87` |
| `data-preparation/3-rename-genes/NEW_Step_5_Replace_gff_Names_with_Dmelanogaster_1_9.py` | `cyp-project-forensics/to_organize/` (same name) | `e974bc122c9f` | `4f93d14a1609` |
| `data-analysis/1-pair-te-with-cyp/repeatOpp.py` | `cyp-project-forensics/pipeline-scripts-output/repeatOpp.py` | `96d5c9eaeb0e` | `3a33ff63f3f7` |
| `data-analysis/1-pair-te-with-cyp/Locate_TE.py` | `cyp-project-forensics/pipeline-scripts-output/Locate_TE.py` | `2ca2029579a3` | `10db4dc9b1f4` |
| `data-analysis/1-pair-te-with-cyp/CleanAnnasse.py` | `cyp-project-forensics/pipeline-scripts-output/CleanAnnasse.py` | `66b306d967a7` | `b4178a1fff98` |
| `data-analysis/1-pair-te-with-cyp/Reg_Gene_Full.txt` | `cyp-project-forensics/pipeline-scripts-output/AnalysisForAll/Reg_Gene_Full.txt` | `df8d7a5f6a12` | unchanged (no header) |
| `data-analysis/2-combine-and-compare/build_tfbs_te_gff.py` | `cyp-project-forensics/analysis-pipeline/build_tfbs_te_gff 1.py` | `21e0453dcc12` | `ef39fd25dd22` |
| `data-analysis/2-combine-and-compare/compare_te_cyp_exposure.py` | `cyp-project-forensics/analysis-pipeline/compare_te_cyp_exposure 1 1.py` | `0cf676ec9f4c` | `fca2c8944460` |
| `data-analysis/2-combine-and-compare/compare_te_cyp_cncc.py` | `cyp-project-forensics/analysis-pipeline/compare_te_cyp_cncc 1 1.py` | `ae5affb95a58` | `fdb9de6cf422` |
| `data-analysis/2-combine-and-compare/compare_te_cyp_xenobiotic.py` | `cyp-project-forensics/analysis-pipeline/compare_te_cyp_xenobiotic 1 1.py` | `5282449d95f6` | `2c4d96bf669c` |

The last four were **renamed** on copy (` 1` / ` 1 1` suffixes dropped): `compare_te_cyp_cncc`
and `compare_te_cyp_xenobiotic` do `import compare_te_cyp_exposure`, which cannot work with a
space in the filename. Apart from the documentation header, contents are unchanged.

The automation's own documentation is copied alongside: `docs/repeat-modeler-automation/README.md`
(unchanged) and `SOURCE.md`. `SOURCE.md` changed in one way only, on 2026-11-24: its `worker.sh:N`
and `rm-manager.sh:N` line references were renumbered to the files as they now are (+34 for
`worker.sh`, whose references are all below its header; +23 for `rm-manager.sh`), so that they
still point at the code they describe.

## Copied data — `datasets/`

| Here | Original | SHA-256 |
|---|---|---|
| `datasets/hog-table/HOG_OG_association_gene_names_without_duplicates_10_31.tsv` | `cyp-project-forensics/to_organize/hog_og/` (same name) | `621473f9336d` |
| `datasets/ananassae-reference/DROSOPHILA_ANANASSAE_final.gff.gz` | Zenodo 18453526 `annotations.tar.gz` → `gffs/` | `c4c45316db82` |
| `datasets/ananassae-reference/Drosophila_ananassae.GCF_017639315.1.rm.fna.out` | `cyp-project-forensics/pipeline-scripts-output/DA_Files/` | `15a67b0c7ce5` |
| `datasets/ananassae-reference/filtered.gff` | `cyp-project-forensics/pipeline-scripts-output/` | `103aae2c7b4a` |
| `datasets/ananassae-reference/GenesAffectedByTEs.txt` | `cyp-project-forensics/pipeline-scripts-output/` | `b094e1a1e22f` |
| `datasets/ananassae-reference/DAnasse_TE_Cyp.txt` | `cyp-project-forensics/pipeline-scripts-output/` | `53d76f762d82` |

---

## Reconstructed: `data-preparation/3-rename-genes/config.py`

`NEW_Step_5_…` and `Step_2_…` both `from config import …`. The lab's `config.py` was never
kept — it is not in the repository, `to_organize/`, or any of its `.7z` archives. This one
supplies what they import:

| Name | Value | Basis |
|---|---|---|
| `STEP5_INPUT` | the `_10_31` HOG table | the only candidate that reproduces the lab's results (below) |
| `DMEL_COLUMN_INDEX` | index of the `DROSOPHILA_MELANOGASTER` column | read from the table header |
| `get_gff_input_path(label)` | `<input>/<label>_final.gff` | Zenodo's filenames |
| `get_gff_output_path(label)` | `<output>/<label>_final_withDmelNames.gff` | the lab's filename for this script's output (lab log; committed example) |
| `clean_field(x)` | strip whitespace, then quotes | quoted multi-gene cells need it; and a whitespace strip explains why Step 2's `_1_9` output lost every line break |
| `STEP2_INPUT`, `STEP2_OUTPUT` | placeholder filenames | Step 2 is reference only |

### Validation

Run on the Zenodo *D. ananassae* annotation with the `_10_31` table:

- **`repeatOpp.py` on the output reproduces the lab's committed `filtered.gff` exactly** —
  all 166 Cyp gene lines (after removing Windows line endings).
- The lab's committed `DA_Files/DROSOPHILA_ANANASSAE_final_withDmelNames.gff` is truncated
  (first 4,329 lines only, ending mid-line). Against those lines, the new output differs only
  on 225 mRNA lines, all explained by the input: the lab's copy of the annotation had tabs
  where Zenodo has spaces inside the description text, so the lab's run kept only the text
  before the first tab. The renaming on those lines (`Parent_New=…`) is identical.

HOG table candidates:

| Table | Result |
|---|---|
| `…_without_duplicates_10_31.tsv` | reproduces the lab exactly — **used** |
| `…_without_duplicates_1_9.tsv` | no line breaks at all; unparseable |
| `…_without_duplicates.tsv` (undated) | many lookups fail; abandoned partway as not matching |

---

## Why `NEW_Step_5_…`, not `ReVamp_Final.py`

The lab had two scripts for the same step, and never recorded which was current.

| | `NEW_Step_5_…` — **used here** | `ReVamp_Final.py` |
|---|---|---|
| Output | `<SPECIES>_final_withDmelNames.gff` | `change_<SPECIES>_final_final.gff` |
| Method | parses attributes; writes `Name=<Dmel>`, keeps `Name_old=` | plain text replacement of IDs/names |
| Multi-gene HOG cells (`"gene-A, gene-B"`) | handled | **every gene in the cell left unrenamed** |
| Short names | exact matches | substring replacement corrupts text (`genome` → `genom,nom,ouibe` in *D. melanogaster*) |
| *D. ananassae* Cyp loci reaching the TE table | 91 | 57 |

What the lab's finished tables (`AnalysisForAll/output/`) were built from, recovered from the
tables themselves: **19 species from ReVamp's `final_final` files** (aldrichi, algonquin,
anomalata, arawakana, arizonae, athabasca, elegans, helvetica, mauritiana, mayaguana, miranda,
pandora, santomea, the four sulfurigaster tables, suzukii, tropicalis), and **7 from a correct
parse** (erecta, mojavensis, pseudoobscura, sechellia, simulans, subpulchrella, willistoni).
In the 19, the 516 Cyp loci ReVamp missed contribute zero TE rows; in the 7, loci of the same
kind contribute 35–178 rows per species. Re-running those 19 species through this pipeline
would add those loci back.

---

## Not included, and why

| File | Reason |
|---|---|
| `ReVamp_Final.py`, `ReVamp_Final_copy_1_10.py`, `make_Dmel_output.py`, `Step5.py`, `Part5_Final.py` | the superseded renaming route |
| `Step_3_Count_CYP_Genes_…`, `Step_4_Extract_CYP_Genes.py` | a side branch (Cyp copy-number tables), not on the path to the comparison |
| `spinContainer.sh`, `runMasker.sh` | the lab's manual RepeatModeler/RepeatMasker runs, replaced by the automation in step 1 |
| `Dmel_HOG_association.tsv`, the INI and gene-list example files | never kept by the lab |
| data (annotations, genomes, HOG table) | too large to copy; the docs say where each comes from |
