#!/usr/bin/env bash
# Stage 2 runner (new, not a lab file, 2026-11-22): renames every species' genes with D. melanogaster
# ortholog names (data preparation step 3) by running the repo's NEW_Step_5 script UNCHANGED.
#
# NEW_Step_5 processes every <SPECIES>_final.gff in GFF_INPUT_DIR, one at a time (about 10 minutes
# each), so this runner parallelises it. It deals the annotations in input/gffs/ out into NB batch
# folders of symlinks: largest file first, round-robin, so the batches are of similar total size.
# It then starts one NEW_Step_5 process per batch, through config.py's environment overrides:
#   STEP5_INPUT=input/hog-table/HOG_OG_association_gene_names_without_duplicates_10_31.tsv
#   GFF_INPUT_DIR=work/batch-<n>   GFF_OUTPUT_DIR=../3-pair-te-with-cyp/input/renamed-gff
# Each process still walks every HOG-table column and prints "Skipping ..." for the species that are
# not in its batch.
# OUTPUT: <SPECIES>_final_withDmelNames.gff in renamed-gff/, and logs/batch-<n>.log per batch. It
#   prints one rc line when all batches are done (rc=1 if any batch failed).
# Usage: bash run-rename.sh [NB]   (default 6; run from anywhere, paths are absolute)
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
