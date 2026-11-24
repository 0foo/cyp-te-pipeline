#!/usr/bin/env bash
# Stage 4b runner (new, not a lab file): runs the three repo compare scripts UNCHANGED, arguments only.
# Needs matplotlib for the plots: run inside nix-shell -p 'python3.withPackages(p: [p.matplotlib])'.
set -euo pipefail
H=/bioinformatics/atallah-repeats/4-combine-and-compare; S=/root/projects/cyp-te-pipeline/data-analysis/2-combine-and-compare; O=$H/output/compare; mkdir -p "$O"; cd "$O"
python3 -B "$S/compare_te_cyp_exposure.py" --self-test
for s in exposure cncc xenobiotic; do
  extra=(); [ $s = xenobiotic ] && extra=(--gene-list "$H/xenobiotic_resistance_cyp_genes.txt")
  echo "== $s"
  python3 -B "$S/compare_te_cyp_$s.py" --config "$H/species_config.ini" "${extra[@]}" \
    --output-csv "$O/$s.species-summary.csv" --per-gene-csv "$O/$s.per-gene.csv" \
    --output-report "$O/$s.report.md" --plot --plot-path "$O/$s.by-species.png"
done
