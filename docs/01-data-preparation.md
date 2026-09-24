# Data preparation

Two files per species come out of this half, and both are needed by data analysis:

| Output | Made by | Used by |
|---|---|---|
| `<Species>.rm.out` — every TE copy in the genome | step 1 | `Locate_TE.py` |
| `<SPECIES>_final_withDmelNames.gff` — genes labelled with *D. melanogaster* ortholog names | step 3 | `repeatOpp.py`, `build_tfbs_te_gff.py` |

The two branches are independent; run them in either order or at the same time.

**Work on local disk.** Downloads, extraction and runs should all happen on a local drive.
Copy results to network storage afterwards if you want them kept there.

---

## Step 1 — repeat library and TE locations

Folder: `data-preparation/1-repeat-library-and-masking/`

**In:** gzipped genome FASTAs (`*.fna.gz`). **Out:** `<sample>.rm.out` per genome.
**Runtime:** 8–26 hours per genome (longer with `LTRSTRUCT=1`).

Per genome, `worker.sh` runs `BuildDatabase` + `RepeatModeler` (builds the species' own repeat
library) and then `RepeatMasker -lib <that library>` (finds every copy in the genome), inside
the `dfam/tetools` container. It is crash-tolerant and runs several genomes in parallel.

```bash
cd data-preparation/1-repeat-library-and-masking
docker pull dfam/tetools:latest
cp rmodeler.conf.example rmodeler.conf
$EDITOR rmodeler.conf        # set IN_DIR, WORK_DIR, OUT_DIR, STATE_DIR, LOG_DIR; keep RUN_MASKER=1
./worker.sh                  # one worker in the foreground first, to see errors
./rm-manager.sh start        # then the real run
./rm-manager.sh status
```

Full reference for the automation (restarts, signals, sizing): the original repository's
`README.md` and `docs/SOURCE.md` in `~/projects/repeat-modeler-automation/`.

The lab originally ran the same two tools by hand (`spinContainer.sh`, `runMasker.sh`); the
automation replaces that and produces the same `.out` format.

---

## Step 2 — the orthogroup (HOG) table — reference only

Folder: `data-preparation/2-orthogroup-table/`

**Do not run these.** They record how the table used in step 3 was made:

- `Step_1_Replace_RNA_identifiers_with_gene_names_Redo.py` turns the orthogroup table's mRNA
  IDs into gene identifiers, reading each species' annotation. Its input
  `Dmel_HOG_association.tsv` was never kept.
- `Step_2_Remove_Duplicate_gene_names_10_31.py` removes duplicate genes within each cell.
  Run with a whitespace-stripping `clean_field` (as the lab's later version did), it strips the
  newline off the last column and **writes every row onto one line** — which is exactly what
  happened to the lab's `_1_9` table.

**Use the table they produced:**

```
~/projects/cyp-project-forensics/to_organize/hog_og/HOG_OG_association_gene_names_without_duplicates_10_31.tsv
```

One row per hierarchical orthogroup (HOG), one column per species. That is `config.py`'s
default. Do not use the `_1_9` table (unparseable) or the undated one (does not reproduce the
lab's results).

---

## Step 3 — rename genes with *D. melanogaster* ortholog names

Folder: `data-preparation/3-rename-genes/`

**In:** the Zenodo annotation for each species + the HOG table.
**Out:** `<SPECIES>_final_withDmelNames.gff`. **Runtime:** about 10 minutes per species.

Why: each species names its genes differently (`LOC6500252`, `gene-G00000000064`…). Giving
every gene the name of its *D. melanogaster* ortholog (`Cyp12e1`) is what lets the Cyp gene
list match across species.

### 3a. Get the annotations

Zenodo record [18453526](https://zenodo.org/records/18453526), 1.76 GB. Zenodo serves it at
about 1–2.5 MB/s whatever you do, so allow 15–25 minutes.

```bash
mkdir -p ~/dl-staging && cd ~/dl-staging
curl -L -C - --retry 10 -o annotations.tar.gz \
  "https://zenodo.org/api/records/18453526/files/annotations.tar.gz/content"
md5sum annotations.tar.gz                   # must be d7cd2d6d0b98b4d51036b05c619c590b
mkdir -p annotations && tar xzf annotations.tar.gz -C annotations
```

That gives `annotations/gffs/<SPECIES>_final.gff.gz` for 301 species.

### 3b. Put the species you want in `input/`

The script processes **every** species whose file is in the input folder, so put in only the
ones you want:

```bash
cd ~/projects/cyp-te-pipeline/data-preparation/3-rename-genes
mkdir -p input
gunzip -c ~/dl-staging/annotations/gffs/DROSOPHILA_ANANASSAE_final.gff.gz \
  > input/DROSOPHILA_ANANASSAE_final.gff
```

Species names are the HOG table's column headers (`DROSOPHILA_ANANASSAE`,
`DROSOPHILA_SULFURIGASTER_BILIMBATA`, …).

### 3c. Run

```bash
python3 NEW_Step_5_Replace_gff_Names_with_Dmelanogaster_1_9.py
```

It prints `on: <SPECIES>` for every column of the HOG table — including the ~300 it skips with
`GFF file not found` — and writes `output/<SPECIES>_final_withDmelNames.gff`.

Paths come from `config.py`; override any of them without editing it:

```bash
STEP5_INPUT=/path/to/hog_table.tsv GFF_INPUT_DIR=/path/to/gffs GFF_OUTPUT_DIR=/path/to/out \
  python3 NEW_Step_5_Replace_gff_Names_with_Dmelanogaster_1_9.py
```

### What the output looks like

Gene lines get the Dmel name first and keep the original as `Name_old`; mRNA lines get
`Parent_New`:

```
…	gene	…	Name=Cyp12e1;ID=gene-G00000000064;Name_old=LOC6500252;…
…	mRNA	…	ID=Dananassae_M00000000004;Parent_New=wek,CG17568;Parent=gene-G00000000004;…
```

Genes with no ortholog, or whose orthogroup has no *D. melanogaster* gene, are left as they
were. *D. melanogaster* itself (and `MUSCA_DOMESTICA`) is skipped — it needs no renaming.

### Checking it worked

Run `repeatOpp.py` (data analysis step 1) on the output. For *D. ananassae* the resulting
`filtered.gff` must be identical to the lab's committed
`~/projects/cyp-project-forensics/pipeline-scripts-output/filtered.gff` (166 lines) once its Windows line
endings are removed:

```bash
cmp <(tr -d '\r' < ~/projects/cyp-project-forensics/pipeline-scripts-output/filtered.gff) filtered.gff
```
