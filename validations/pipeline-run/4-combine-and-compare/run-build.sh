#!/usr/bin/env bash
# Stage 4a runner (new, not a lab file, 2026-11-23): builds one combined GFF3 per species (data analysis
# step 2a) by running the repo's build_tfbs_te_gff.py UNCHANGED, with arguments only. It covers every
# species in exposure.tsv, which is all 301 stage 3 species.
#   --te-hits         input/te-tables/<LABEL>_GenesAffectedByTE.txt   (from stage 3, run-pairing.sh)
#   --annotation-gff  the renamed GFF (stage 2); D. melanogaster uses its own un-renamed GFF
#   --genome-fasta    work/<LABEL>/genome.fna: the unmasked genome, decompressed here because the
#                     script indexes plain FASTA by byte offset. It is deleted afterwards.
#   --fimo-via-docker --fimo-docker-image memesuite-root:local --jaspar-cache-dir work/jaspar
#       FIMO runs in memesuite/memesuite rebuilt to run as root (fimo-image/Dockerfile). The official
#       image runs as user meme, which cannot read the script's mode-700 temp dir, so FIMO failed and the
#       script carried on without binding sites. A run whose log says 'FIMO run failed' counts as FAIL.
#   --output          output/combined-gff/<LABEL>.combined.gff3. The script's defaults also write
#                     .transposable_elements.gff3 and bgzip/tabix copies.
# A failed species' output is deleted. Species whose output already exists are skipped, so the runner
# can be re-run after an interruption.
# OUTPUT: one OK/FAIL/SKIP line per species (with time and feature count). Logs go to logs/run-build.log
#   and to logs/build/<LABEL>.log, one per species.
# Usage: bash run-build.sh [JOBS]   (default 6)
set -uo pipefail
JOBS=${1:-6}
export H=/bioinformatics/atallah-repeats/4-combine-and-compare
export S=/root/projects/cyp-te-pipeline/data-analysis/2-combine-and-compare
export G=/mnt/wd/bioinformatics/genomes/2_unmasked_datasets
mkdir -p "$H/output/combined-gff" "$H/logs/build" "$H/work"
one() {
    label=$1; rm_sample=$2
    out=$H/output/combined-gff/$label.combined.gff3
    [ -s "$out" ] && { echo "$(date '+%T') SKIP $label (exists)"; return; }
    w=$H/work/$label; mkdir -p "$w"
    if [ "$label" = DROSOPHILA_MELANOGASTER ]; then
        gff=/bioinformatics/atallah-repeats/2-rename-genes/input/gffs/DROSOPHILA_MELANOGASTER_final.gff
    else
        gff=/bioinformatics/atallah-repeats/3-pair-te-with-cyp/input/renamed-gff/${label}_final_withDmelNames.gff
    fi
    t0=$(date +%s)
    zcat "$G/$rm_sample.rm.fna.gz" > "$w/genome.fna" \
    && (cd "$w" && python3 -B "$S/build_tfbs_te_gff.py" --te-hits "$H/input/te-tables/${label}_GenesAffectedByTE.txt" \
            --annotation-gff "$gff" --genome-fasta "$w/genome.fna" --fimo-via-docker --fimo-docker-image memesuite-root:local \
            --jaspar-cache-dir "$H/work/jaspar" --output "$out") > "$H/logs/build/$label.log" 2>&1 \
    && ! grep -q 'FIMO run failed' "$H/logs/build/$label.log" \
    && echo "$(date '+%T') OK   $label $(( $(date +%s) - t0 ))s $(grep -vc '^#' "$out") features" \
    || { rm -f "$out"; echo "$(date '+%T') FAIL $label (see logs/build/$label.log)"; }
    rm -f "$w/genome.fna" "$w/genome.fna.fai"
}
export -f one
awk -F'\t' 'NR>1{print $1"\t"$2}' "$H/exposure.tsv" \
    | xargs -P "$JOBS" -d '\n' -I{} bash -c 'IFS=$'"'"'\t'"'"' read -r l s <<< "{}"; one "$l" "$s"'
echo "$(date '+%F %T') done"
