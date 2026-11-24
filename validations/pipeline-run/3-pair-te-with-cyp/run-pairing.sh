#!/usr/bin/env bash
# Stage 3 runner (new, not a lab file, 2026-11-22): pairs Cyp genes with nearby repeats (data analysis
# step 1) for every species, by running the three lab scripts repeatOpp.py -> Locate_TE.py ->
# CleanAnnasse.py.
#
# The lab scripts hardcode Windows paths at module scope, and CLAUDE.md forbids editing them. So for
# each row of pairs.tsv with status RUN, this runner makes runs/<GFF_LABEL>/ and writes path-patched
# COPIES of the three scripts there, replacing only the hardcoded-path lines with sed. The lab copies
# in the repo are never edited. It then runs the copies in that directory, because they write fixed
# filenames into the current directory. Up to JOBS species run in parallel (xargs -P).
#   repeatOpp.py  line 46 GFF = <renamed GFF> (D. melanogaster: its own un-renamed GFF);
#                 line 48 rgFile = Reg_Gene_Full.txt
#   Locate_TE.py  line 49 cypGene = filtered.gff; line 52 TEs = input/rm-out/<rm_sample>.rm.out
#   CleanAnnasse.py line 43 aliasToFind = Reg_Gene_Full.txt   (toClean is already relative)
# These are the line numbers of the lab scripts as they now stand, with their 2026-11-24 documentation
# header. In the header-free lab originals the same lines were 4/6, 2/5 and 3, which is what the full
# 301-species run used. Re-check them if a header changes length (docs/02-data-analysis.md).
# OUTPUT: runs/<LABEL>/{filtered.gff, GenesAffectedByTEs.txt, DAnasse_TE_Cyp.txt, run.log}. The final
#   table is copied to ../4-combine-and-compare/input/te-tables/<GFF_LABEL>_GenesAffectedByTE.txt.
#   One OK/FAIL line is printed per species.
# Needs pandas (imported, unused, by repeatOpp.py); run inside:
#   nix-shell -p 'python3.withPackages(p: [p.pandas])' --run 'bash run-pairing.sh [JOBS]'   (default 8)
set -euo pipefail
JOBS=${1:-8}
A=/bioinformatics/atallah-repeats/3-pair-te-with-cyp
export A REPO=/root/projects/cyp-te-pipeline/data-analysis/1-pair-te-with-cyp
export OUT=/bioinformatics/atallah-repeats/4-combine-and-compare/input/te-tables
mkdir -p "$OUT" "$A/runs" "$A/logs"
one() {
    rm_sample=$1; label=$2
    w=$A/runs/$label; mkdir -p "$w"; cd "$w"
    if [ "$label" = DROSOPHILA_MELANOGASTER ]; then
        gff=/bioinformatics/atallah-repeats/2-rename-genes/input/gffs/DROSOPHILA_MELANOGASTER_final.gff
    else
        gff=$A/input/renamed-gff/${label}_final_withDmelNames.gff
    fi
    sed -e "46s|.*|GFF = \"$gff\"|" -e "48s|.*|rgFile =\"$REPO/Reg_Gene_Full.txt\"|" "$REPO/repeatOpp.py" > repeatOpp.py
    sed -e "49s|.*|cypGene = \"filtered.gff\"|" -e "52s|.*|TEs = '$A/input/rm-out/$rm_sample.rm.out'|" "$REPO/Locate_TE.py" > Locate_TE.py
    sed -e "43s|.*|aliasToFind = '$REPO/Reg_Gene_Full.txt'|" "$REPO/CleanAnnasse.py" > CleanAnnasse.py
    { python3 repeatOpp.py && python3 Locate_TE.py && python3 CleanAnnasse.py; } > run.log 2>&1 \
        && cp DAnasse_TE_Cyp.txt "$OUT/${label}_GenesAffectedByTE.txt" \
        && echo "$(date '+%T') OK   $label $(wc -l < filtered.gff) cyp-lines $(wc -l < DAnasse_TE_Cyp.txt) rows" \
        || echo "$(date '+%T') FAIL $label (see runs/$label/run.log)"
}
export -f one
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="status") c=i; next} $c=="RUN"{print $1"\t"$2}' "$A/pairs.tsv" \
    | xargs -P "$JOBS" -d '\n' -I{} bash -c 'IFS=$'"'"'\t'"'"' read -r s l <<< "{}"; one "$s" "$l"'
echo "$(date '+%F %T') done"
