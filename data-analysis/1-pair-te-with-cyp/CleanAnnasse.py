# ===== DOCUMENTATION HEADER (added 2026-11-24; not part of the lab original) =====
# CleanAnnasse.py: data analysis, step 1c. Replaces the gene-attribute column of the TE table with a
# single Cyp name.
#
# PIPELINE POSITION
#   The last of repeatOpp.py -> Locate_TE.py -> CleanAnnasse.py, in data-analysis/1-pair-te-with-cyp/.
#   Its output is the "TE hits" table that build_tfbs_te_gff.py (data analysis step 2a) reads.
#   The name comes from the lab's first test species, D. ananassae ("Annasse"); the script works
#   unchanged for any species.
#
# HOW IT WORKS
#   1. Reads the Cyp name list (aliasToFind, Reg_Gene_Full.txt), strips each line, and skips blank
#      lines and lines starting with "D" (the 5 non-melanogaster IDs). That leaves 92 names, in file order.
#   2. Splits each line of Locate_TE.py's output on tabs. Lines with fewer than 2 fields are kept
#      as they are. Otherwise column 2 (the gene's full GFF attribute text) is replaced by the first
#      list name that occurs in it as a substring. If no name occurs, column 2 is left unchanged.
#   3. Writes every line, tab-joined, to DAnasse_TE_Cyp.txt.
#
# WHY
#   Downstream scripts key everything on this name. build_tfbs_te_gff.py looks it up in the annotation,
#   tags each TE "Within range of <name>", and the compare_te_cyp_*.py scripts count TEs per name.
#
# INPUTS (hardcoded at module scope; see docs/02-data-analysis.md)
#   toClean      the "toClean = ..." line below: GenesAffectedByTEs.txt (already a relative path)
#   aliasToFind  the "aliasToFind = ..." line below: Reg_Gene_Full.txt from this folder (must be edited)
# OUTPUT
#   DAnasse_TE_Cyp.txt, a fixed name in the current directory: the same rows with column 2 = one Cyp name.
#   Rename it to D_<species>GenesAffectedByTE.txt. run-pairing.sh copies it to
#   <LABEL>_GenesAffectedByTE.txt instead. D. ananassae check: identical to
#   datasets/ananassae-reference/DAnasse_TE_Cyp.txt (900 rows, 465 distinct) apart from CRLF.
#
# KNOWN DEFECTS (kept as-is, because this is an archival copy)
#   - Substring matching, first match in list order. When names nest (Cyp4g1 inside Cyp4g15), or
#     when a locus carries several Dmel names (Name=Cyp313a5,Cyp313a2,...), the label is whichever
#     list entry comes first. That is not necessarily the gene's own name, and paralogous loci can
#     end up with the same label.
#   - Duplicate rows from repeatOpp.py pass straight through.
# Run in this project by validations/pipeline-run/3-pair-te-with-cyp/run-pairing.sh, which writes a
# path-patched copy per species (sed on the aliasToFind = line only) and leaves this file untouched.
# ===== END DOCUMENTATION HEADER =====
toClean = 'GenesAffectedByTEs.txt'

aliasToFind = 'C:/Users/User/Documents/BioAtallah/RepeatMasker/RepeatOpp/Dmel_Only/Dmel_Data/Reg_Gene_Full.txt'

with open(toClean, mode = 'r',encoding = 'utf-8') as file, open(aliasToFind, mode = 'r', encoding = 'utf-8') as names, open('DAnasse_TE_Cyp.txt', mode = 'w', encoding = 'utf-8') as out: 
    nameList = []
    container = []

    for name in names:
        nm = name.strip()
        if not nm:
            continue
        if nm.startswith("D"):
            continue
        nameList.append(nm)
    #print(nameList)


    for line in file:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 2:
            container.append(fields)
            continue

        for name in nameList:
            if name in fields[1]:
                fields[1] = name
                break
        container.append(fields)

    for ele in container:
        out.write("\t".join(ele) + "\n")