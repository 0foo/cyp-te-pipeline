# ===== DOCUMENTATION HEADER (added 2026-11-24; not part of the lab original) =====
# repeatOpp.py: data analysis, step 1a. Selects the Cyp gene records from one species' renamed
# annotation.
#
# PIPELINE POSITION
#   The first of three scripts in data-analysis/1-pair-te-with-cyp/, run once per species in a
#   scratch directory: repeatOpp.py -> Locate_TE.py -> CleanAnnasse.py. Input comes from data
#   preparation step 3 (NEW_Step_5). Its output, filtered.gff, is the input of Locate_TE.py.
#
# HOW IT WORKS
#   1. Reads the Cyp name list (rgFile, Reg_Gene_Full.txt, one name per line). Lines starting with
#      "D" are skipped: these are the 5 non-melanogaster IDs such as Dvir\GJ21722. That leaves 92
#      D. melanogaster Cyp names.
#   2. Reads the GFF one line at a time, skips "#" lines, and keeps each line whose column 3 is
#      "gene" and whose column 9 (the attribute text) contains any list name as a plain substring.
#   3. Writes the kept lines, unchanged, to filtered.gff.
#
# THEORY
#   Gene IDs differ between species (LOC6500252, gene-G00000000064, ...). NEW_Step_5 has written the
#   D. melanogaster ortholog's name into each gene's Name= attribute, using the HOG orthogroup table.
#   That is what lets a single list of D. melanogaster Cyp names select the Cyp genes of any species.
#
# INPUTS (hardcoded at module scope; edit them before every run, see docs/02-data-analysis.md)
#   GFF     the "GFF = ..." line below: <SPECIES>_final_withDmelNames.gff (GFF3, tab-separated)
#   rgFile  the "rgFile = ..." line below: Reg_Gene_Full.txt from this folder
# OUTPUT
#   filtered.gff, a fixed name written to the current directory: the selected gene lines, verbatim.
#   Expected result for D. ananassae: 166 lines (91 distinct), identical to
#   datasets/ananassae-reference/filtered.gff once its CRLF line endings are removed.
#
# KNOWN DEFECTS (kept as-is, because this is an archival copy)
#   - No `break` in the name loop. A gene line whose attributes contain several list names is
#     appended once per match, so filtered.gff has duplicate lines and each copy later gets its own
#     TE rows. build_tfbs_te_gff.py removes the duplicates again.
#   - Substring matching. Cyp4g1 also matches Cyp4g15, and a name can match anywhere in the
#     attributes (Name_old=, a description, and so on).
#   - Genes that the renamed annotation names by symbol (the Halloween genes phm, dib, sad, shd) do not
#     contain the "Cyp306a1"-style numbers the list uses, so they are missed except in D. melanogaster.
#   - `import pandas as pd` is never used, but pandas (2.x) must still be installed.
# Run in this project by validations/pipeline-run/3-pair-te-with-cyp/run-pairing.sh, which writes a
# path-patched copy per species (sed on the GFF = and rgFile = lines) and leaves this file untouched.
# ===== END DOCUMENTATION HEADER =====
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