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