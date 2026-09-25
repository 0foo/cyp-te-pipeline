# cyp-te-pipeline

The happy path through the Atallah lab's Cyp / transposable-element pipeline: which
transposable elements (TEs) sit in or near cytochrome P450 (Cyp) genes, species by species,
and whether that differs between species with high and low xenobiotic exposure.

Every script here is a **byte-for-byte copy** of a lab original (checksums in
[`docs/03-provenance.md`](docs/03-provenance.md)). The one exception is
`data-preparation/3-rename-genes/config.py`, a reconstruction of a lost lab file, labelled as
such and validated against the lab's own results.

The repository is self-contained **for the validated *D. ananassae* path**: the orthogroup table
and a complete *D. ananassae* reference set are in [`datasets/`](datasets/), and no download is
needed to run and check that path. Going beyond it needs two separate downloads from Zenodo
record [18453526](https://zenodo.org/records/18453526): the annotation archive
`annotations.tar.gz` (1.76 GB) for species other than *D. ananassae*, and the genome archive
`genomes.tar.gz` (19.6 GB) for the repeat-masking step and for `build_tfbs_te_gff.py`. Single
genomes can be taken from NCBI instead — *D. ananassae* is GCF_017639315.1.

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

| Folder | Contents | Known issues |
|---|---|---|
| [`data-preparation/1-repeat-library-and-masking/`](data-preparation/1-repeat-library-and-masking/) | `worker.sh`, `rm-manager.sh`, `rmodeler.conf.example` — genome → repeat library → TE locations (`.rm.out`) | — |
| [`data-preparation/2-orthogroup-table/`](data-preparation/2-orthogroup-table/) | `Step_1_…`, `Step_2_…` — how the HOG table was built. **Reference only**: they cannot be re-run, and the table they made is the input to step 3 | — |
| [`data-preparation/3-rename-genes/`](data-preparation/3-rename-genes/) | `NEW_Step_5_…` + reconstructed `config.py` — relabels each species' genes with *D. melanogaster* ortholog names | `config.py` is a reconstruction, not a lab original |
| [`data-analysis/1-pair-te-with-cyp/`](data-analysis/1-pair-te-with-cyp/) | `repeatOpp.py`, `Locate_TE.py`, `CleanAnnasse.py`, `Reg_Gene_Full.txt` — Cyp genes × nearby TEs | **[no repeat filtering, edge-straddling TEs discarded, duplicate rows, substring matching](#known-issues)** |
| [`data-analysis/2-combine-and-compare/`](data-analysis/2-combine-and-compare/) | `build_tfbs_te_gff.py`, `compare_te_cyp_{exposure,cncc,xenobiotic}.py` — combined annotation, then the cross-species statistics | **[passes unfiltered repeats into the statistics](#known-issues)** |
| [`datasets/`](datasets/) | the HOG table, and the *D. ananassae* inputs and lab outputs used to check the pipeline | — |
| [`docs/`](docs/) | how to run each step, and where every file came from | — |
| [`NOTES.md`](NOTES.md) | open gaps: what is still missing or unverified | — |

## Start here

1. [`docs/01-data-preparation.md`](docs/01-data-preparation.md) — getting a species from raw
   downloads to a TE table and a renamed annotation.
2. [`docs/02-data-analysis.md`](docs/02-data-analysis.md) — from those two files to the
   cross-species comparison.
3. [`docs/03-provenance.md`](docs/03-provenance.md) — what was copied from where, what was
   reconstructed, how it was validated, and why this renaming script and not the other one.

## Why this is the "correct" path

The lab had two gene-renaming scripts. This pipeline uses
`NEW_Step_5_Replace_gff_Names_with_Dmelanogaster_1_9.py`, not `ReVamp_Final.py`,
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

## Known issues

### The TE hits are mostly not transposable elements

`Locate_TE.py` applies no filtering beyond position. It keeps a repeat if it is on the same
chromosome as the gene and lies entirely within the gene extended by 3,000 bp — there is no
length floor, no Smith-Waterman score threshold, no divergence limit and no repeat-class
exclusion. `build_tfbs_te_gff.py` de-duplicates rows but does not filter by class, length or
divergence either, so everything reaches the cross-species statistics.

Measured on the validated *D. ananassae* output — `GenesAffectedByTEs.txt`, 900 rows:

| Repeat class | Rows | Share of rows | Median length |
|---|---|---|---|
| `Simple_repeat` | 502 | 55.8% | 39 bp |
| `Low_complexity` | 118 | 13.1% | 46 bp |
| `Unknown` (unclassified) | 136 | 15.1% | 161 bp |
| Classified TE families | 144 | 16.0% | — |

The median hit is **47 bp**, the shortest **6 bp**. Half the rows are under 50 bp and 78% are
under 100 bp. The five shortest hits are `(ACAC)n` and `(GCCCTC)n` microsatellites of 6–8 bp.
Smith-Waterman scores run from **11** with a median of **19**; alignment divergence has a median
of 16.4%.

**68.9% of rows (620) are `Simple_repeat` or `Low_complexity`, which are not transposable
elements** — they are microsatellite and low-complexity sequence that RepeatMasker reports
alongside genuine repeat families. Those 620 rows total just **27,945 bp**.

Reducing the table to what is actually a classified transposable element:

| | Value |
|---|---|
| Rows in the table | 900 |
| Distinct rows | 521 |
| Distinct repeat placements (chromosome, start, end, family) | 457 |
| Distinct Cyp gene labels | 77 |
| Rows with a classified TE family | 144 |
| Distinct TE placements behind those rows | **63** |
| Cyp gene labels with any classified TE nearby | **29 of 77** |

So a table of 900 rows rests on **63 distinct transposable-element placements** affecting
**29 Cyp gene labels**. Adding a length floor on top of the class filter:

| Filter | Rows retained |
|---|---|
| none (current behaviour) | 900 |
| classified TE family | 144 |
| classified TE family, ≥ 80 bp | 100 |
| classified TE family, ≥ 100 bp | 87 |

**Consequence.** The raw row count should not be reported as a count of TEs near Cyp genes; as
produced it is closer to a microsatellite-density measurement. Any result taken from this
pipeline needs either the class breakdown alongside it or an explicit statement of the filter
applied.

Filtering is not to be added to the scripts themselves — they are archival copies. It belongs in
a post-processing step outside this repository, or stated as a caveat on the results.

### Transposons straddling the window edge are discarded whole

`Locate_TE.py` requires **both ends** of a repeat to fall inside the search window:

```python
start = int(cypFields[3]) - 3000          # gene start, minus flank
stop  = int(cypFields[4]) + 3000          # gene end,   plus flank
if cypFields[0] == eleFields[4] and (int(eleFields[6]) <= stop and int(eleFields[5]) >= start):
```

The window is the gene body plus 3,000 bp on each side, so both flanks are searched. But a
repeat that begins inside the window and continues past its edge fails `te_end <= stop` and is
dropped **entirely** — not truncated, not recorded as a partial overlap, simply absent. The same
happens at the left edge via `te_begin >= start`. The conventional test is overlap
(`te_begin <= stop and te_end >= start`), which keeps it.

Nothing about the Cyp gene is lost — the gene was already selected by `repeatOpp.py`. What is
lost is the *transposon*, so the gene appears to have fewer elements near it than it does.

**This biases the table towards short repeats.** A repeat of length *L* placed in a 3,000 bp
flank survives only if it ends before the flank does, so its chance of being captured is
`max(0, 1 - L/3000)`:

| Repeat length | Chance of capture in a flank |
|---|---|
| 50 bp | 98.3% |
| 100 bp | 96.7% |
| 500 bp | 83.3% |
| 1,000 bp | 66.7% |
| 1,500 bp | 50.0% |
| 2,000 bp | 33.3% |
| 2,500 bp | 16.7% |
| **≥ 3,000 bp** | **0% — cannot be captured in a flank at all** |

In this genome's `.rm.out` (297,070 repeat annotations, 142,835 of them classified TE families),
**4,096 classified TEs — 2.9% — are ≥ 3,000 bp and can therefore never be found in a flank**;
13.9% are ≥ 1,000 bp and so lose at least a third of their chance. The effect falls hardest on
exactly the families of interest: the median `LINE/L2` in this genome is 1,355 bp (55% ceiling)
and the median `DNA/TcMar` is 923 bp (69% ceiling), whereas a 39 bp microsatellite is captured
essentially always.

**Measured loss on *D. ananassae*.** Re-running the same pairing with an overlap test instead of
a containment test:

| | Containment (current) | Overlap | Discarded |
|---|---|---|---|
| Rows | 900 | 930 | 30 |
| Median row length | 47 bp | — | **1,150 bp** |
| Rows with a classified TE family | 144 | 156 | 12 |
| Median length of those | 137 bp | — | **1,961 bp** |
| Distinct classified-TE placements | 63 | 69 | **6** |

Only 30 rows are lost — 3.3% — but they are 24× longer at the median than what is kept, and they
carry **6 of the 69 distinct transposable-element placements, 8.7% of the real signal**. Of the
30, 16 had at least 100 bp inside the window and 7 had at least 1,000 bp inside.

**Transposons discarded despite lying substantially inside the window:**

| Family | Class | Length | Inside window | Edge crossed | Cyp gene |
|---|---|---|---|---|---|
| `ltr-1_family-11` | `LTR/Pao` | 6,339 bp | 1,555 bp (24.5%) | right | `Cyp313a5,Cyp313a2,Cyp313a3,Cyp313a1` |
| `rnd-4_family-78` | `DNA/hAT-Ac` | 2,129 bp | 1,756 bp (**82.5%**) | left | `Cyp6g1` |
| `rnd-1_family-520` | `RC/Helentron` | 1,793 bp | 14 bp (0.8%) | right | `Cyp309a1,Cyp309a2` |
| `rnd-1_family-25` | `LINE/CR1` | 1,150 bp | 228 bp (19.8%) | left | `Cyp9b1,Cyp9b2` |
| `rnd-1_family-520` | `RC/Helentron` | 478 bp | 161 bp (33.7%) | right | `Cyp4c3` |
| `rnd-4_family-107` | `LINE/R1` | 176 bp | 144 bp (**81.8%**) | right | `Cyp28c1` |

Coordinates for the first two, so they can be checked directly:

- `NC_057927.1:4270106-4276444` against window `4263690-4271660` — a 6,339 bp `LTR/Pao` element
  whose first 1,555 bp sit inside the window of the `Cyp313a` cluster. Discarded.
- `NC_057929.1:5461814-5463942` against window `5462187-5470756` — a 2,129 bp `DNA/hAT-Ac`
  element, 82.5% of it inside `Cyp6g1`'s window, overlapping the gene's upstream flank.
  Discarded. `Cyp6g1` is the locus most often cited in *Drosophila* insecticide resistance, and
  the element is missing from the table purely because its first 373 bp lie outside the cut-off.

**Consequence.** Because the loss is length-dependent rather than random, it cannot be treated as
noise. Long insertions — full-length LTR elements and LINEs, the ones most likely to carry
regulatory sequence — are the ones systematically absent, and they are absent in a way no count
of the output reveals. Combined with the class composition above, the table over-represents
microsatellites twice over: they dominate what is kept, and they are the only class that the
containment rule almost never rejects.

### Reproducing the lab's outputs does not validate the measurement

The *D. ananassae* check confirms the code still behaves as it did originally: `filtered.gff`
(166 Cyp gene rows), `GenesAffectedByTEs.txt` (900 rows, 521 distinct) and `DAnasse_TE_Cyp.txt`
(900 rows, 465 distinct) are reproduced byte-for-byte. That is a test of faithfulness, not of validity — it would pass
equally well if the underlying measurement were meaningless. The issue above is not detectable
by this check, nor is the boundary defect above, nor either quirk below.

### Duplicate rows inflate apparent hit counts

`repeatOpp.py` has no `break` in its name-matching loop, so a gene whose attributes contain
several Cyp names is emitted once per match, and each copy accumulates its own TE rows. The
*D. ananassae* `DAnasse_TE_Cyp.txt` is 900 rows but only 465 distinct; the 900 rows resolve to
457 distinct repeat placements. Row counts are therefore not hit counts.

### Substring matching mislabels nested gene names

Cyp name matching is a plain substring test, so `Cyp4g1` also matches `Cyp4g15`. `CleanAnnasse.py`
labels each row with the *first* list entry found in the attribute text, which is not necessarily
the gene's own name wherever names nest.

### Unfiltered comparisons may track assembly provenance rather than biology

Simple-repeat and low-complexity content varies with assembly quality and sequencing platform,
and the available genomes are a mix of long-read and short-read assemblies. Because the
comparison scripts receive unfiltered repeat hits, a high- versus low-exposure difference could
reflect how the genomes were sequenced and assembled rather than any biological signal. This
should be resolved before any cross-species result is presented — alongside the missing exposure
assignments, which are not recorded anywhere in the lab's files.
