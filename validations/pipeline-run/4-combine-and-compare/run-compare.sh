#!/usr/bin/env bash
# Stage 4b runner (new, not a lab file, 2026-11-23): runs the three repo comparison scripts UNCHANGED,
# with arguments only. That is data analysis step 2b: the high- versus low-exposure test of TEs near
# Cyp genes.
#   1. compare_te_cyp_exposure.py --self-test   (the statistics cross-checks; the script stops on failure)
#   2. For s in exposure, cncc and xenobiotic: compare_te_cyp_<s>.py --config species_config.ini
#      (from make-exposure.py). The xenobiotic script also gets --gene-list
#      xenobiotic_resistance_cyp_genes.txt (10 genes).
# OUTPUTS, in 4-combine-and-compare/output/compare/: <s>.species-summary.csv, <s>.per-gene.csv,
#   <s>.report.md and <s>.by-species.png. sheets/make-sheet.py and make-hits-sheet.py read them.
# Caveat: the scripts count every repeat row as a "TE", simple repeats included, and pool genes across
#   species. See the headers of the compare scripts and README.md, "Known issues".
# Needs matplotlib for the plots: run inside nix-shell -p 'python3.withPackages(p: [p.matplotlib])',
#   then bash run-compare.sh.
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
