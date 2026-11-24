# ===== DOCUMENTATION HEADER (added 2026-11-24; not part of the lab original) =====
# Locate_TE.py: data analysis, step 1b. Pairs each selected Cyp gene with the repeats near it.
#
# PIPELINE POSITION
#   The second of repeatOpp.py -> Locate_TE.py -> CleanAnnasse.py, in data-analysis/1-pair-te-with-cyp/.
#   It reads repeatOpp.py's filtered.gff and the species' RepeatMasker table, which comes from data
#   preparation step 1 (worker.sh: a RepeatModeler de novo library, then RepeatMasker -lib).
#
# HOW IT WORKS
#   1. Loads every non-blank line of the RepeatMasker .out file into memory. The 3 header lines are
#      loaded too, but never match a chromosome name.
#   2. For each gene line of filtered.gff, builds a window running from (gene start - 3000) to
#      (gene end + 3000): the gene body plus 3 kb of flank on each side, regardless of strand.
#   3. Tests every .out row against the window. It keeps a row when the query sequence (.out column 5)
#      equals the gene's seqid and the repeat lies entirely inside the window, i.e.
#      repeat end <= window end AND repeat begin >= window start.
#   4. Writes one line per (gene line, repeat) pair: seqid TAB gene attributes (GFF column 9) TAB the
#      complete .out row, with its original spacing.
#   The search is a plain nested loop (genes x repeats, with no index). It takes seconds to minutes.
#
# THEORY
#   A TE inserted in a gene, or in the few kb around it, can change that gene's expression or
#   structure; a classic case is the Accord insertion upstream of Cyp6g1 that gives DDT resistance.
#   The 3 kb flank is the lab's operational definition of "near".
#
# INPUTS (hardcoded at module scope; edit them before every run, see docs/02-data-analysis.md)
#   cypGene  the "cypGene = ..." line below: filtered.gff from repeatOpp.py
#   TEs      the "TEs = ..." line below: <species>.rm.out, a RepeatMasker .out table
#            (whitespace-separated: SW score, %div, %del, %ins, query, begin, end, (left), strand,
#            repeat, class/family, repeat positions, ID)
# OUTPUT
#   GenesAffectedByTEs.txt, a fixed name in the current directory, tab-separated as described above.
#   Expected result for D. ananassae: identical to datasets/ananassae-reference/GenesAffectedByTEs.txt
#   (900 rows, 521 distinct) apart from CRLF line endings.
#
# KNOWN DEFECTS AND LIMITATIONS (kept as-is, because this is an archival copy; see README.md, "Known issues")
#   - No filtering beyond position: no repeat-class, length, score or divergence filter. In D. ananassae,
#     68.9% of rows are Simple_repeat or Low_complexity (microsatellites, which are not TEs), and the
#     median hit is 47 bp.
#   - The containment test drops any repeat that straddles a window edge, rather than keeping the
#     overlap. That biases the table against long elements: a repeat of 3 kb or more can never be
#     found in a flank. In D. ananassae, 30 rows are lost, including 6 of 69 distinct TE placements
#     (for example a 2,129 bp hAT element that lies 82.5% inside Cyp6g1's window).
#   - The duplicate gene lines from repeatOpp.py are carried through, each with its own copy of the TE rows.
# Run in this project by validations/pipeline-run/3-pair-te-with-cyp/run-pairing.sh, which writes a
# path-patched copy per species (sed on the cypGene = and TEs = lines) and leaves this file untouched.
# ===== END DOCUMENTATION HEADER =====
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

