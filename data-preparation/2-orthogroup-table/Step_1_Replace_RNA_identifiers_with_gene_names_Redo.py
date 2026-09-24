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
    
 





