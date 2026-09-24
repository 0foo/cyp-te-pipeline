import pandas as pd

#Gff file
GFF = "C:/Users/User/Documents/BioAtallah/RepeatMasker/RepeatOpp/DA_Files/DROSOPHILA_ANANASSAE_final_withDmelNames.gff"

rgFile ="C:/Users/User/Documents/BioAtallah/RepeatMasker/RepeatOpp/Dmel_Only/Dmel_Data/Reg_Gene_Full.txt"

stripGFF = []
regGene = []


with open(rgFile, mode = "r", encoding="utf-8") as t:
    for gene in t:
        if gene.startswith('D'):
           continue
        regGene.append(gene.rstrip('\n'))
        #print(gene)

#print(regGene)
        

with open(GFF, mode = "r", encoding="utf-8") as s:

    for line in s:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")

        for Rgene in regGene:
            #print(Rgene)
            if fields[2] == "gene" and Rgene in fields[8]:
                stripGFF.append(line)
    #print(stripGFF)


with open("filtered.gff", "w", encoding="utf-8") as out:
    for row in stripGFF:
        out.write(row)