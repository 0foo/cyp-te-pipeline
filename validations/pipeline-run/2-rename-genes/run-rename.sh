#!/usr/bin/env bash
# Stage 2 runner (new, not a lab file). Splits the annotations into NB batches, each its own
# GFF_INPUT_DIR of symlinks to input/gffs/, and runs the unmodified NEW_Step_5 script once per
# batch, in parallel. Output goes to ../3-pair-te-with-cyp/input/renamed-gff/.
# Usage: bash run-rename.sh [NB]   (default 6)
set -euo pipefail
NB=${1:-6}
S=/bioinformatics/atallah-repeats/2-rename-genes
REPO=/root/projects/cyp-te-pipeline/data-preparation/3-rename-genes
OUT=/bioinformatics/atallah-repeats/3-pair-te-with-cyp/input/renamed-gff
HOG=$S/input/hog-table/HOG_OG_association_gene_names_without_duplicates_10_31.tsv
mkdir -p "$S/work" "$S/logs" "$OUT"
rm -rf "$S"/work/batch-*
# deal the files out largest-first, round-robin, so batches have similar total size
i=0
while read -r f; do
    b=$(( i % NB )); mkdir -p "$S/work/batch-$b"; ln -s "$f" "$S/work/batch-$b/"; i=$((i+1))
done < <(ls -S "$S"/input/gffs/*_final.gff)
cd "$REPO"
pids=()
for d in "$S"/work/batch-*; do
    b=$(basename "$d")
    STEP5_INPUT=$HOG GFF_INPUT_DIR=$d GFF_OUTPUT_DIR=$OUT \
        python3 -u NEW_Step_5_Replace_gff_Names_with_Dmelanogaster_1_9.py > "$S/logs/$b.log" 2>&1 &
    pids+=($!)
done
rc=0; for p in "${pids[@]}"; do wait "$p" || rc=1; done
echo "$(date '+%F %T') all batches finished rc=$rc"
