#cyp genes that are affected by TEs
cypGene = "C:/Users/User/Documents/BioAtallah/RepeatMasker/RepeatOpp/filtered.gff"

#All possible TEs
TEs = 'C:/Users/User/Documents/BioAtallah/RepeatMasker/RepeatOpp/DA_Files/Drosophila_ananassae.GCF_017639315.1.rm.fna.out'

#print(TEs)


with open(cypGene, mode="r", encoding="utf-8") as cyp,open(TEs, mode="r", encoding="utf-8") as te, open("GenesAffectedByTEs.txt", mode="w", encoding="utf-8") as out:
    # preload TE file lines so we can iterate multiple times safely
    te_lines = [line.rstrip('\n') for line in te if line.strip()]

    for gene in cyp:
        gene = gene.strip()

        if not gene:
            continue

        cypFields = gene.split("\t")
        #print(cypFields)
        if len(cypFields) == 0:
            continue

        for ele in te_lines:
            eleFields = ele.split()
            # ensure expected column exists before indexing
            if len(eleFields) <= 4:
                continue
                
            start = int(cypFields[3]) - 3000
            stop = int(cypFields[4]) + 3000
            


            #print(eleFields[5])
            if cypFields[0] == eleFields[4] and (int(eleFields[6])<= stop and int(eleFields[5]) >= start):
               
                line = f"{cypFields[0]}\t{cypFields[8]}\t{ele}\n"
                out.write(line) 

