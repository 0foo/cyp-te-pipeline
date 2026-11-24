# ===== DOCUMENTATION HEADER (added 2026-11-24; not part of the lab original) =====
# Step_1_Replace_RNA_identifiers_with_gene_names_Redo.py: data preparation, step 2 (orthogroup
# table). REFERENCE ONLY: it cannot be re-run, because its input Dmel_HOG_association.tsv was never kept.
#
# WHAT IT DID
#   It converted the orthogroup table from transcript (mRNA) IDs to gene identifiers. That is the
#   first of two steps that produced the HOG table which data preparation step 3 (NEW_Step_5) uses to
#   rename genes. The table actually used is datasets/hog-table/..._without_duplicates_10_31.tsv.
#
# HOW IT WORKS
#   1. Reads HOG_OG_association_1_9.tsv: one row per hierarchical orthogroup (HOG), with columns
#      HOG, OG, gene-tree parent clade, and then one column per species listing that species' mRNA IDs
#      separated by ", ". The header row holds the species labels (DROSOPHILA_ANANASSAE, ...).
#   2. Reads Dmel_HOG_association.tsv and maps each HOG to a D. melanogaster common gene name (its column 5).
#   3. For each species column (skipping the first 3 columns, MUSCA_DOMESTICA and DROSOPHILA_MELANOGASTER),
#      reads annotations/gff_fixed/<SPECIES>_final.gff and builds a map from mRNA ID to gene identifier,
#      using the mRNA's Name= value, or its Parent= value when Name= is missing or holds a RefSeq "XM" accession.
#   4. Rewrites every cell: each mRNA ID becomes its gene identifier, or "<mRNA>_NOT_FOUND". The
#      D. melanogaster column is filled from step 2, and the OG and clade columns are copied.
#   5. Writes Output/HOG_OG_association_gene_names_Redo_1_9.tsv: same layout, tab-separated, cells
#      ", "-joined, with no newline after the last row.
#
# THEORY
#   Orthogroups group genes that descend from one ancestral gene. Isoforms of one gene therefore
#   become repeated gene names within a cell, and Step_2 removes them. Gene-level names are needed
#   because the renaming step edits gene features. The orthology tool is not recorded; the HOG/OG
#   columns look like OrthoFinder output, from the "Dhakad_et_al_2025_Data" release.
#
# INPUTS / OUTPUTS: relative paths hardcoded at the top (Dhakad_et_al_2025_Data/..., Output/...).
# LIMITATIONS: subprocess and os are imported but unused. The file handles are opened at module scope.
#   The date suffixes (_1_9 here, _10_31 on the table in use) do not line up; see docs/03-provenance.md.
# ===== END DOCUMENTATION HEADER =====
import subprocess
import os
import re

HOG_OG_association_lines = open("Dhakad_et_al_2025_Data/HOG_OG_association_1_9.tsv").readlines()
Dmel_HOG_association_lines = open("Dhakad_et_al_2025_Data/Dmel_HOG_association.tsv").readlines()
output = open("Output/HOG_OG_association_gene_names_Redo_1_9.tsv", "w")

#The first line of the HOG_OG_association file consists of labels, including species names
labels = HOG_OG_association_lines[0].rstrip("\n").split("\t")

#Create a dictionary of HOGs from the input (mRNA) file
HOG_mRNAs = {}
#Create another dictionary of HOGs with gene names
HOG_genes = {}
#Create a dictionary for Dmel genes
Dmel_HOG_genes = {}

for i in range(1,len(HOG_OG_association_lines)):
    fields = (HOG_OG_association_lines[i]).rstrip("\n").split("\t")
    #Create a dictionary where the key is the HOG and the value is a list of the associated mRNAs in each species
    HOG_mRNAs[fields[0]] = []
    for j in range(1,len(fields)):
        HOG_mRNAs[fields[0]].append(fields[j])

for i in range(1,len(Dmel_HOG_association_lines)):
    fields = (Dmel_HOG_association_lines[i]).rstrip("\n").split("\t")
    #Create a dictionary where the key is the HOG and the value is the Dmel common name
    Dmel_HOG_genes[fields[0]] = fields[4]    

#This loop goes through the HOG_OG_association file column by column (working on one species at a time)
for i in range(1, len(labels)):
    mRNA_to_gene_dictionary = {}
    #Print the column
    print(labels[i])
    if ((i>2) and (labels[i] != "MUSCA_DOMESTICA") and (labels[i] != "DROSOPHILA_MELANOGASTER")): #There is no gff file associated with "MUSCA_DOMESTICA". There are also no OGs associated with this species in the table (all blank)        
        gff_file_name = "Dhakad_et_al_2025_Data/annotations/gff_fixed/" + labels[i] + "_final.gff"
        gff_file = open(gff_file_name)
        for line in gff_file:
            line = line.rstrip("\n")
            fields = line.split("\t")
            if(len(fields)>7):
                if (fields[2] == "mRNA"):
                    description = fields[8]
                    #print(description)
                    #match = re.search(r"ID=([a-zA-Z0-9_-\.]*)", description)
                    match = re.search(r"ID=([^;]*)", description)
                    mRNA = match.group(1)
                    #print(mRNA)
                    gene = ""
                    if ("Name=" in description):
                        #match = re.search(r"Name=([a-zA-Z0-9_-\.]*)", description)
                        match = re.search(r"Name=([^;]*)", description)
                        gene = match.group(1)
                        if ("XM" in gene):
                            #match = re.search(r"Parent=([a-zA-Z0-9_-\.]*)", description)
                            match = re.search(r"Parent=([^;]*)", description)
                            gene = match.group(1)
                        else:
                            gene = match.group(1)
                    else:
                        #match = re.search(r"Parent=([a-zA-Z0-9_-\.]*)", description)
                        match = re.search(r"Parent=([^;]*)", description)
                        gene = match.group(1)
                    #print(gene)
                    mRNA_to_gene_dictionary[mRNA] = gene
        gff_file.close()        
    for HOG in HOG_mRNAs.keys():
        if(i<2):
            HOG_genes[HOG] = []
        if (i<3):
            HOG_genes[HOG].append(HOG_mRNAs[HOG][i-1])
        else:
            if (labels[i] == "DROSOPHILA_MELANOGASTER"):
                #HOG_genes[HOG].append("")
                HOG_genes[HOG].append(Dmel_HOG_genes.get(HOG,""))
            else:    
                mRNAs = HOG_mRNAs[HOG][i-1].split(", ")
                #print(mRNAs)
                for j in range(0, len(mRNAs)):
                    if (j==0):
                        HOG_genes[HOG].append("")
                    else:
                        HOG_genes[HOG][i-1] = HOG_genes[HOG][i-1] + ", "
                    if (mRNAs[j] != ""):
                        HOG_genes[HOG][i-1] = HOG_genes[HOG][i-1] + mRNA_to_gene_dictionary.get(mRNAs[j], mRNAs[j] + "_NOT_FOUND")
 
#Write labels to output file
output.write(HOG_OG_association_lines[0].rstrip("\n"))

for HOG in HOG_genes.keys():
    output.write("\n")
    output.write(HOG)
    for gene in HOG_genes[HOG]:
        output.write("\t")
        output.write(gene)

output.close()        
    
 





