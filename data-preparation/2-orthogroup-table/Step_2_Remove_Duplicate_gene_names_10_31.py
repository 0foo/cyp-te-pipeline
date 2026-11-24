# ===== DOCUMENTATION HEADER (added 2026-11-24; not part of the lab original) =====
# Step_2_Remove_Duplicate_gene_names_10_31.py: data preparation, step 2 (orthogroup table), second
# half. REFERENCE ONLY: do not run it. Use the table it produced,
# datasets/hog-table/HOG_OG_association_gene_names_without_duplicates_10_31.tsv.
#
# WHAT IT DOES
#   It removes repeated gene names within each species cell of the HOG table written by Step_1.
#   Several isoforms (mRNAs) of one gene map to the same gene name, so after Step_1 a cell can list
#   one gene several times.
#
# HOW IT WORKS
#   1. Reads STEP2_INPUT, which is imported from config.py. The header line is copied unchanged.
#   2. For each data row, splits the line on tabs. The first 3 columns (HOG, OG, gene-tree parent clade)
#      are copied. Each later cell is split on ",", every part goes through config.clean_field() (strip
#      whitespace, then quotes), and the parts are de-duplicated in first-seen order and joined with "," (no space).
#   3. Writes STEP2_OUTPUT and prints "duplicate found: <gene>" for every duplicate removed.
#
# INPUTS / OUTPUTS
#   STEP2_INPUT, STEP2_OUTPUT and clean_field come from config.py, which is a reconstruction (the lab's
#   was lost), so here they are placeholder filenames. Both files are tab-separated HOG tables: one row
#   per HOG, one column per species.
#
# KNOWN DEFECT (kept as-is, because this is an archival copy)
#   The row is split without stripping its trailing newline, so the newline belongs to the last cell.
#   A whitespace-stripping clean_field removes it, and every row is then written onto one line. That
#   is exactly what happened to the lab's unusable _1_9 table. The _10_31 table still has its line
#   breaks, so it must have been written with a clean_field that left the newline alone.
#   subprocess, os and re are imported but unused.
# ===== END DOCUMENTATION HEADER =====
import subprocess
import os
import re
from config import STEP2_INPUT, STEP2_OUTPUT, clean_field

HOG_OG_association_gene_names_file = open(STEP2_INPUT, "r")
output = open(STEP2_OUTPUT, "w")

HOG_OG_association_gene_names_lines = HOG_OG_association_gene_names_file.readlines()

labels = HOG_OG_association_gene_names_lines[0].rstrip("\n").split("\t")

label_line = HOG_OG_association_gene_names_lines[0]
#print(label_line)

for line_index, line in enumerate(HOG_OG_association_gene_names_lines):
    
    if (line_index == 0):
        # Write header line
        output.write(label_line)
    
    else: 
        values = line.split("\t")
        for i in range(0, len(labels)):
            if (i == 0):
                HOG_Value = values[i]
                output.write(HOG_Value)
            if (i == 1):
                OG_Value = values[i]
                output.write("\t")
                output.write(OG_Value)
            if (i == 2):      
                Gene_Parent_Tree_Clade_Value = values[i]
                output.write("\t")
                output.write(Gene_Parent_Tree_Clade_Value)
            if (i > 2):
                output.write("\t")
                GeneUsed = []  # Use list to preserve order
                GeneUsed_set = set()  # Use set for fast lookup
                Genes = values[i]
                gene_names = [clean_field(x) for x in Genes.split(",")]
                for z in gene_names:
                    if (z not in GeneUsed_set):
                        GeneUsed.append(z)
                        GeneUsed_set.add(z)
                    else:
                        print("duplicate found:", z)
                GeneUsed_counter = 0
                for gene in GeneUsed:
                    if (GeneUsed_counter > 0):
                        output.write(",")
                    output.write(gene)
                    GeneUsed_counter = GeneUsed_counter + 1


HOG_OG_association_gene_names_file.close()
output.close()