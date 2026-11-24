# cyp-genes-from-hog-table

**New (not lab data, 2026-11-24).** Every cytochrome P450 (Cyp) gene that the orthogroup table
`../hog-table/HOG_OG_association_gene_names_without_duplicates_10_31.tsv` assigns to each species. This is
the pipeline's Cyp gene set **before renaming**: the genes stage 2 (`NEW_Step_5`) is supposed to rename to
*D. melanogaster* names, and stage 3 then searches for.

Rebuild with `python3 -B extract_cyp_genes.py` (standard library only; the method and caveats are in its header).

## How a gene counts as a Cyp

A HOG row (one hierarchical orthogroup) counts as a Cyp orthogroup when its `DROSOPHILA_MELANOGASTER` cell
holds a P450 name. That means either `Cyp<number><letters>…` or one of the Halloween P450s that FlyBase names
by symbol: phm, dib, sad, shd, spo, spok. 65 rows qualify, and none of them mixes in non-P450 names.
Every species gene in those rows is listed here. Membership comes from the orthology inference alone, not
from sequence or domain evidence; the inference tool was never recorded, though the naming looks like OrthoFinder.

## Files (tab-separated, header row)

| File | Rows | Contents |
|---|---|---|
| `cyp_genes.tsv` | 25,177 genes | `species` · `gene_id_in_hog_table` (exactly as in the table) · `id_type` · `dmel_cyp` (the D. melanogaster Cyp, or Cyp set, of the orthogroup) · `n_dmel_cyps_in_hog` · `halloween` (symbol and Cyp number, if a Halloween gene) · `hog` · `og` · `gene_tree_parent_clade` · `genes_in_cell` (that species' gene count in this orthogroup) |
| `cyp_counts_per_species.tsv` | 305 species | Cyp genes; distinct D. melanogaster Cyps covered; Cyp HOGs with ≥1 gene; stable Halloween genes present (of 4: phm, dib, sad, shd); number of IDs of each `id_type` |
| `cyp_copy_number_matrix.tsv` | 305 × 65 | species × D. melanogaster Cyp orthogroup → number of that species' genes in it |

`id_type` records how the table wrote the identifier, and this matters downstream:

| id_type | Example | Genes | Matches in the species' GFF |
|---|---|---|---|
| `gene_id` | `gene-G00000020176` | 8,698 | `ID=`; the only form NEW_Step_5 looks up when a gene has an ID |
| `ncbi_locus_name` | `LOC111079515` | 14,275 | `Name=` / `gene_name=` of genes projected from an NCBI reference |
| `annotation_gene_id` | `AMIOTA_COM_G0010083` | 1,801 | `gene_id=` / `Name=` |
| `gene_symbol` | `Cyp6a21`, `phm` | 387 | mostly the D. melanogaster column itself (85) |
| `refseq_transcript` | `NM_001364801.1` | 10 | transcript accessions |
| `other` | `CG33964` | 6 | FlyBase CG numbers |

Rows with any `id_type` other than `gene_id` are the ones stage 2 fails to rename. That is about 65% of
Cyp genes in the median species, and it is why the scripts later find so few Cyps
(`validations/README.md`, and the "Cyp-TE validation" sheet).

## Summary (305 species columns)

| Measure | Value |
|---|---|
| Cyp genes per species | median 81 (range 0–173) |
| D. melanogaster Cyps covered | median 78 of 85 |
| Stable Halloween genes present | median 4 of 4 |
| Species with no Cyp genes | *Musca domestica* (only one species) |

Expected P450 count per Drosophila genome is 74–94 (Good et al. 2014, *GBE* 6:1118,
[doi:10.1093/gbe/evu083](https://doi.org/10.1093/gbe/evu083)), so the table's orthology captures about the
full set.

## Caveats

- **Multi-Cyp orthogroups.** In 15 orthogroups the D. melanogaster cell lists several Cyps (e.g. `Cyp9f3,Cyp9f2`).
  Their species genes cannot be assigned to one of them, so `dmel_cyp` keeps the whole set and the matrix counts it as one column.
  As a result, "D. melanogaster Cyps covered" can exceed the number of genes.
- **Missing D. melanogaster Cyps.** Cyp6a9, Cyp6a14, Cyp6a16 and Cyp6a20 are in no orthogroup, so no species has them here.
- **Extra species columns.** The table has 305 species columns; the pipeline ran 301 species (stage 3).
  The others are americana, paucipunta, pseudoobscura_pseudoobscura and *Musca*, which have no annotation.
