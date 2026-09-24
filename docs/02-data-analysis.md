# Data analysis

**In, per species:** `<SPECIES>_final_withDmelNames.gff` and `<Species>.rm.out` from data
preparation, plus the genome FASTA. **Out:** a CSV, a markdown report and a plot comparing TE
burden near Cyp genes between high- and low-exposure species.

The copied scripts keep the lab's **hardcoded paths**. Nothing here edits them for you: before
each run, open the script and set the lines listed below. Work in a scratch directory per
species — the step 1 scripts write fixed filenames into the current directory.

---

## Step 1 — pair Cyp genes with nearby TEs

Folder: `data-analysis/1-pair-te-with-cyp/`. Three scripts, run in order, per species.
**Runtime:** seconds to a couple of minutes.

| # | Script | Set these lines first | Reads | Writes (current directory) |
|---|---|---|---|---|
| 1 | `repeatOpp.py` | 4 `GFF` → the `_withDmelNames.gff`; 6 `rgFile` → `Reg_Gene_Full.txt` | renamed GFF + Cyp list | `filtered.gff` |
| 2 | `Locate_TE.py` | 2 `cypGene` → `filtered.gff`; 5 `TEs` → the species' `.rm.out` | `filtered.gff` + `.out` | `GenesAffectedByTEs.txt` |
| 3 | `CleanAnnasse.py` | 1 `toClean` → `GenesAffectedByTEs.txt`; 3 `aliasToFind` → `Reg_Gene_Full.txt` | the above + Cyp list | `DAnasse_TE_Cyp.txt` |

```bash
REPO=/path/to/cyp-te-pipeline
mkdir -p runs/ananassae && cd runs/ananassae      # any scratch directory
S=$REPO/data-analysis/1-pair-te-with-cyp
python3 $S/repeatOpp.py && python3 $S/Locate_TE.py && python3 $S/CleanAnnasse.py
mv DAnasse_TE_Cyp.txt D_ananassaeGenesAffectedByTE.txt
```

The output file is always called `DAnasse_TE_Cyp.txt` (it was written for *D. ananassae*);
**rename it** to `D_<species>GenesAffectedByTE.txt`.

**What each does:**

1. **`repeatOpp.py`** keeps the `gene` lines whose attributes contain any name in
   `Reg_Gene_Full.txt` — 92 *D. melanogaster* Cyp gene names. (The file has 97 entries; the 5
   starting with `D`, such as `Dvir\GJ21722`, are non-*melanogaster* IDs and both scripts skip
   them.) This is where the renaming pays off:
   without `Name=Cyp12e1` the gene would not match.
2. **`Locate_TE.py`** keeps every TE on the same chromosome that lies **entirely within** the
   gene extended by **3,000 bp** on each side. One output row per (gene, TE) pair: chromosome,
   the gene's attributes, then the full RepeatMasker line.
3. **`CleanAnnasse.py`** replaces the attribute text in column 2 with the first Cyp name from
   the list found in it.

**Behaviour to know about** (lab code, kept as-is):

- **Duplicate rows.** `repeatOpp.py` has no `break`, so a gene line matching several Cyp names
  is written once per match, and each copy gets its own set of TE rows. The committed
  ananassae table has 900 rows, 465 unique. `build_tfbs_te_gff.py` removes the duplicates.
- **Substring matching.** `Cyp4g1` also matches `Cyp4g15`. Column 2's label is the first list
  entry found, which may not be the gene's own name when names nest.

**Check:** everything needed is in `datasets/ananassae-reference/`. Point `Locate_TE.py` at
its `filtered.gff` and `Drosophila_ananassae.GCF_017639315.1.rm.fna.out`, run steps 2–3, and
the outputs match its `GenesAffectedByTEs.txt` and `DAnasse_TE_Cyp.txt` exactly (apart from
CRLF line endings).

---

## Step 2 — combine and compare

Folder: `data-analysis/2-combine-and-compare/`. These scripts take command-line arguments; no
path editing needed. The files were renamed on copy (the lab's had ` 1` / ` 1 1` suffixes) so
that the two restricted comparison scripts can `import compare_te_cyp_exposure`.

### 2a. One combined annotation per species — `build_tfbs_te_gff.py`

Merges the step 1 table, the renamed annotation and the genome into one GFF3: Cyp genes,
their transcripts/exons/introns, the TEs, and predicted transcription-factor binding sites
(JASPAR insect motifs plus the CncC:Maf-S motif, scanned with FIMO).

```bash
S=$REPO/data-analysis/2-combine-and-compare
python3 $S/build_tfbs_te_gff.py \
  --te-hits D_ananassaeGenesAffectedByTE.txt \
  --annotation-gff DROSOPHILA_ANANASSAE_final_withDmelNames.gff \
  --genome-fasta Drosophila_ananassae.GCF_017639315.1.fna \
  --output combined_ananassae.gff3
```

- Needs `fimo` (MEME Suite) on `PATH`, or `--fimo-via-docker`, or `--fimo-path`. Add
  `--skip-tfbs` to leave out binding sites — but then `compare_te_cyp_cncc.py` has nothing to
  work with.
- `--sequence-source ncbi --ncbi-email you@…` fetches only the needed sequence from NCBI
  instead of reading a local genome.
- `python3 build_tfbs_te_gff.py --help` documents every flag; the settings can also live in an
  INI file (`--config`, `--species`).

**Runtime:** minutes per species.

### 2b. Cross-species comparison

All three scripts read the same INI file, one section per species:

```ini
[ananassae]
gff      = combined_ananassae.gff3
exposure = high

[santomea]
gff      = combined_santomea.gff3
exposure = low
```

Paths are relative to the INI file. **Which species are high- and low-exposure is not recorded
anywhere in the lab's files** (only *D. santomea* is noted as low); you have to supply it.

```bash
python3 $S/compare_te_cyp_exposure.py   --config species.ini                       # every Cyp gene
python3 $S/compare_te_cyp_cncc.py       --config species.ini                       # only Cyp genes near a CncC:Maf-S site
python3 $S/compare_te_cyp_xenobiotic.py --config species.ini --gene-list genes.txt # a curated resistance list
```

Each writes a summary CSV, a markdown report and a PNG plot (`--no-plot` to skip;
`--per-gene-csv` for per-gene rows). `--self-test` runs each script's built-in checks.

`genes.txt` for the xenobiotic run is one gene name per line; `#` starts a comment. The lab's
example list (`xenobiotic_resistance_cyp_genes.example.txt`) was never kept.

---

## *D. melanogaster*

It needs no renaming (its genes already have their own names), so use its Zenodo
`DROSOPHILA_MELANOGASTER_final.gff` directly as the annotation in both steps.
