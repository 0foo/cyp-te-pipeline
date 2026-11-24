#!/usr/bin/env bash
# Stage 3 runner (new, not a lab file). For every row of pairs.tsv with status RUN, makes
# runs/<GFF_LABEL>/, writes path-patched copies of the three lab scripts there (the lab copies in
# the repo are never edited), runs repeatOpp -> Locate_TE -> CleanAnnasse, and copies the result to
# ../4-combine-and-compare/input/te-tables/<GFF_LABEL>_GenesAffectedByTE.txt.
# Needs pandas (imported, unused, by repeatOpp.py); run inside:
#   nix-shell -p 'python3.withPackages(p: [p.pandas])' --run 'bash run-pairing.sh [JOBS]'
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
    sed -e "4s|.*|GFF = \"$gff\"|" -e "6s|.*|rgFile =\"$REPO/Reg_Gene_Full.txt\"|" "$REPO/repeatOpp.py" > repeatOpp.py
    sed -e "2s|.*|cypGene = \"filtered.gff\"|" -e "5s|.*|TEs = '$A/input/rm-out/$rm_sample.rm.out'|" "$REPO/Locate_TE.py" > Locate_TE.py
    sed -e "3s|.*|aliasToFind = '$REPO/Reg_Gene_Full.txt'|" "$REPO/CleanAnnasse.py" > CleanAnnasse.py
    { python3 repeatOpp.py && python3 Locate_TE.py && python3 CleanAnnasse.py; } > run.log 2>&1 \
        && cp DAnasse_TE_Cyp.txt "$OUT/${label}_GenesAffectedByTE.txt" \
        && echo "$(date '+%T') OK   $label $(wc -l < filtered.gff) cyp-lines $(wc -l < DAnasse_TE_Cyp.txt) rows" \
        || echo "$(date '+%T') FAIL $label (see runs/$label/run.log)"
}
export -f one
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="status") c=i; next} $c=="RUN"{print $1"\t"$2}' "$A/pairs.tsv" \
    | xargs -P "$JOBS" -d '\n' -I{} bash -c 'IFS=$'"'"'\t'"'"' read -r s l <<< "{}"; one "$s" "$l"'
echo "$(date '+%F %T') done"
