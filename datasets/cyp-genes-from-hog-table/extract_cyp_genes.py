#!/usr/bin/env python3
# NEW (not a lab file, 2026-11-24): extracts every cytochrome P450 (Cyp) gene from the orthogroup (HOG) table.
#
# What it does
#   The HOG table (../hog-table/HOG_OG_association_gene_names_without_duplicates_10_31.tsv) has one row per
#   hierarchical orthogroup and one column per species. Its cells list each species' genes in that orthogroup.
#   The D. melanogaster cell gives the FlyBase gene name(s). A HOG row is a Cyp orthogroup when its D. melanogaster
#   cell holds a P450 name, i.e. Cyp<number><letters>[number] (Cyp6g1, Cyp12d1-p...) or one of the six Halloween
#   P450s FlyBase names by symbol (phm, dib, sad, shd, spo, spok). All 65 such rows contain only P450 names. Every
#   species gene listed in those rows is written out. This is the pipeline's Cyp gene set BEFORE renaming (stage 2):
#   the genes NEW_Step_5 is supposed to rename to D. melanogaster names. Identifiers are copied exactly as the table
#   has them (gene IDs like gene-G00000020176, or gene names like LOC111079515 / <SPECIES>_G0005554).
#
# Theory and caveats
#   - Orthology is inferred by the table's source tool, whose name was never recorded (the N1.HOG / OG naming looks
#     like OrthoFinder). A gene is a Cyp here only by orthogroup membership, not by sequence or domain evidence.
#   - In 15 rows the D. melanogaster cell lists several Cyps (e.g. "Cyp9b1,Cyp9b2"). The species genes in such a row
#     cannot be assigned to one of them, so each is recorded with the whole name set. The matrix counts that set as
#     one column, so no gene is counted twice.
#   - 4 D. melanogaster P450s (Cyp6a9, Cyp6a14, Cyp6a16, Cyp6a20) are in no HOG row, so no species has them here.
#
# Inputs: the HOG table (path below, relative to this file). Standard library only.
# Outputs (next to this script, tab-separated, UTF-8, "\n" line endings):
#   cyp_genes.tsv              one row per (species, gene): species, gene_id_in_hog_table, id_type, dmel_cyp,
#                              n_dmel_cyps_in_hog, halloween, hog, og, gene_tree_parent_clade, genes_in_cell
#   cyp_counts_per_species.tsv one row per species: Cyp genes, distinct D. melanogaster Cyps covered, Cyp HOGs with
#                              >=1 gene, the 4 stable Halloween genes present (of 4), and how the IDs are written
#   cyp_copy_number_matrix.tsv species x D. melanogaster Cyp (or Cyp set) -> number of genes in that orthogroup
# Usage: python3 -B extract_cyp_genes.py
import csv, os, re, collections

HERE = os.path.dirname(os.path.abspath(__file__))
HOG = os.path.join(HERE, "..", "hog-table", "HOG_OG_association_gene_names_without_duplicates_10_31.tsv")
P450 = re.compile(r"^Cyp\d+[a-z]+\d*", re.I)
HALLOWEEN = {"phm": "Cyp306a1", "dib": "Cyp302a1", "sad": "Cyp315a1", "shd": "Cyp314a1", "spo": "Cyp307a1", "spok": "Cyp307a2"}
STABLE_HALLOWEEN = ["phm", "dib", "sad", "shd"]   # single-copy in all 12 genomes, Good et al. 2014 Table 1

def id_type(g):
    if g.startswith(("gene-", "rna-", "Gene-", "RNA-")): return "gene_id"         # matches the GFF ID= attribute
    if g.startswith("LOC"): return "ncbi_locus_name"                                   # GFF Name=/gene_name=
    if re.match(r"^(NM|XM|NR|XR)_\d", g): return "refseq_transcript"
    if re.match(r"^[A-Z0-9_]+_G\d+$", g): return "annotation_gene_id"                  # GFF gene_id= / Name=
    if P450.match(g) or g in HALLOWEEN: return "gene_symbol"                                # a FlyBase-style name
    return "other"

csv.field_size_limit(10**9)
rd = csv.reader(open(HOG, newline=""), delimiter="\t")
hdr = [h.strip() for h in next(rd)]
dm = hdr.index("DROSOPHILA_MELANOGASTER")
species = hdr[3:]
genes, labels = [], []
for row in rd:
    row = [c.strip() for c in row]
    names = [n.strip() for n in row[dm].split(",") if n.strip()]
    if not names or not all(P450.match(n) or n in HALLOWEEN for n in names): continue
    label = ",".join(names)
    if label not in labels: labels.append(label)
    for i, sp in enumerate(species, start=3):
        cell = [g.strip() for g in row[i].split(",") if g.strip()] if i < len(row) else []
        for g in cell:
            genes.append(dict(species=sp, gene_id_in_hog_table=g, id_type=id_type(g), dmel_cyp=label,
                              n_dmel_cyps_in_hog=len(names),
                              halloween=",".join(f"{n} ({HALLOWEEN[n]})" for n in names if n in HALLOWEEN),
                              hog=row[0], og=row[1], gene_tree_parent_clade=row[2], genes_in_cell=len(cell)))

def write(name, cols, rows):
    with open(os.path.join(HERE, name), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)

write("cyp_genes.tsv", list(genes[0]), genes)

by_sp = collections.defaultdict(list)
for g in genes: by_sp[g["species"]].append(g)
counts = []
for sp in species:
    gs = by_sp.get(sp, [])
    covered = {n for g in gs for n in g["dmel_cyp"].split(",")}
    ty = collections.Counter(g["id_type"] for g in gs)
    counts.append(dict(species=sp, cyp_genes=len(gs), dmel_cyps_covered=len(covered), cyp_hogs_with_genes=len({g["hog"] for g in gs}),
                       stable_halloween_present=sum(h in covered for h in STABLE_HALLOWEEN),
                       ids_gene_id=ty["gene_id"], ids_ncbi_locus_name=ty["ncbi_locus_name"],
                       ids_annotation_gene_id=ty["annotation_gene_id"], ids_refseq_transcript=ty["refseq_transcript"], ids_gene_symbol=ty["gene_symbol"], ids_other=ty["other"]))
write("cyp_counts_per_species.tsv", list(counts[0]), counts)

mat = []
for sp in species:
    c = collections.Counter(g["dmel_cyp"] for g in by_sp.get(sp, []))
    mat.append({"species": sp, **{l: c.get(l, 0) for l in labels}})
write("cyp_copy_number_matrix.tsv", ["species"] + labels, mat)
print(f"{len(genes)} Cyp genes, {len(species)} species, {len(labels)} D. melanogaster Cyp orthogroup labels")
