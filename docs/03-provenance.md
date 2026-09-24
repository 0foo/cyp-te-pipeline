# Provenance

Where every file in this folder came from, what was reconstructed, and the evidence behind
choosing this path. Assembled 2026-09-24.

The lab's code and data are treated as evidence: nothing in `~/projects/cyp-project` or
`~/projects/repeat-modeler-automation` was changed to build this. The forensic write-up these
findings come from is `~/projects/cyp-project/docs/deep/06-final-final-gff.md`.

---

## Copied files

All byte-for-byte copies; the first 12 hex digits of each SHA-256 match the original.

| Here | Original | SHA-256 |
|---|---|---|
| `data-preparation/1-repeat-library-and-masking/worker.sh` | `repeat-modeler-automation/worker.sh` | `bbbe404c8c80` |
| `data-preparation/1-repeat-library-and-masking/rm-manager.sh` | `repeat-modeler-automation/rm-manager.sh` | `b61435b5e58b` |
| `data-preparation/1-repeat-library-and-masking/rmodeler.conf.example` | `repeat-modeler-automation/rmodeler.conf.example` | `8d551841c8e0` |
| `data-preparation/2-orthogroup-table/Step_1_Replace_RNA_identifiers_with_gene_names_Redo.py` | `cyp-project/to_organize/` (same name) | `74d53b16c85f` |
| `data-preparation/2-orthogroup-table/Step_2_Remove_Duplicate_gene_names_10_31.py` | `cyp-project/to_organize/` (same name) | `1d06b9a8d914` |
| `data-preparation/3-rename-genes/NEW_Step_5_Replace_gff_Names_with_Dmelanogaster_1_9.py` | `cyp-project/to_organize/` (same name) | `e974bc122c9f` |
| `data-analysis/1-pair-te-with-cyp/repeatOpp.py` | `cyp-project/pipeline-scripts-output/repeatOpp.py` | `96d5c9eaeb0e` |
| `data-analysis/1-pair-te-with-cyp/Locate_TE.py` | `cyp-project/pipeline-scripts-output/Locate_TE.py` | `2ca2029579a3` |
| `data-analysis/1-pair-te-with-cyp/CleanAnnasse.py` | `cyp-project/pipeline-scripts-output/CleanAnnasse.py` | `66b306d967a7` |
| `data-analysis/1-pair-te-with-cyp/Reg_Gene_Full.txt` | `cyp-project/pipeline-scripts-output/AnalysisForAll/Reg_Gene_Full.txt` | `df8d7a5f6a12` |
| `data-analysis/2-combine-and-compare/build_tfbs_te_gff.py` | `cyp-project/analysis-pipeline/build_tfbs_te_gff 1.py` | `21e0453dcc12` |
| `data-analysis/2-combine-and-compare/compare_te_cyp_exposure.py` | `cyp-project/analysis-pipeline/compare_te_cyp_exposure 1 1.py` | `0cf676ec9f4c` |
| `data-analysis/2-combine-and-compare/compare_te_cyp_cncc.py` | `cyp-project/analysis-pipeline/compare_te_cyp_cncc 1 1.py` | `ae5affb95a58` |
| `data-analysis/2-combine-and-compare/compare_te_cyp_xenobiotic.py` | `cyp-project/analysis-pipeline/compare_te_cyp_xenobiotic 1 1.py` | `5282449d95f6` |

The last four were **renamed** on copy (` 1` / ` 1 1` suffixes dropped): `compare_te_cyp_cncc`
and `compare_te_cyp_xenobiotic` do `import compare_te_cyp_exposure`, which cannot work with a
space in the filename. Contents are unchanged.

To re-verify:

```bash
cd ~/projects && sha256sum cyp-te-pipeline/data-analysis/1-pair-te-with-cyp/repeatOpp.py \
  cyp-project/pipeline-scripts-output/repeatOpp.py
```

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

## Why Ayush's renaming, not `ReVamp_Final.py`

The lab had two scripts for the same step, and never recorded which was current.

| | `NEW_Step_5_…` (Ayush) — **used here** | `ReVamp_Final.py` (Duy) |
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
