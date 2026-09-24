#!/usr/bin/env python3
"""
compare_te_cyp_xenobiotic.py

Compares TE concentration ACROSS THE 5 SPECIES for a specific, literature-
curated list of Cyp genes known to be involved in xenobiotic/insecticide
resistance (e.g. Cyp6g1 - the classic case where an Accord TE insertion in
its own promoter drives DDT/neonicotinoid resistance) - rather than every
annotated Cyp gene (compare_te_cyp_exposure.py) or every CncC-proximal Cyp
gene (compare_te_cyp_cncc.py).

REQUIRES: compare_te_cyp_exposure.py in the same directory - this script
imports its already-validated GFF3 parsing and statistics functions
rather than re-implementing/duplicating them.

GENE LIST
---------
The target gene list is a plain-text file, NOT hardcoded in this script -
see xenobiotic_resistance_cyp_genes.example.txt for the format and a
starter list with literature citations. That starter list is a STARTING
POINT compiled from well-known Drosophila insecticide-resistance
literature, not an authoritative or exhaustive one - review/edit it
against your own literature search before treating results as final.
Copy the .example file, rename it, and edit freely (plain text, no code
changes needed).

Matching a target symbol against a species' Cyp genes is case-insensitive
and exact (not partial/regex) - if a species' annotation splits a gene
into paralogs (e.g. Cyp12d1 appearing as Cyp12d1-d/Cyp12d1-p), list each
variant explicitly in the gene-list file or it won't match.

INPUT
-----
Same combined GFF3 files and --config file as compare_te_cyp_exposure.py.
Only needs `gene` and TE-Description records - no TF_binding_site
dependency, so (unlike compare_te_cyp_cncc.py) this comparison works for
every species that has usable gene/TE data, regardless of whether FIMO
motif scanning was ever run on it.

OUTPUT
------
Same shape as the other two scripts: per-species CSV, markdown report
(with the same pooled Fisher's exact / Mann-Whitney U tests, now run only
on the target gene list), and an optional bar plot. Also reports, per
species, which target genes were NOT found at all (diagnostic for naming
mismatches across species).

USAGE
-----
    python compare_te_cyp_xenobiotic.py --config te_cyp_species_config.ini \\
        --gene-list xenobiotic_resistance_cyp_genes.txt

    python compare_te_cyp_xenobiotic.py --self-test
"""

import argparse
import csv
import sys
from pathlib import Path

try:
    import compare_te_cyp_exposure as base
except ImportError:
    sys.exit(
        "[error] could not import compare_te_cyp_exposure.py - this script needs it in the "
        "same directory (it reuses its validated GFF3 parsing and statistics functions)."
    )


# ==========================================================================
# 1. Target gene list parsing

def parse_gene_list(path):
    """Reads a plain-text target-gene-list file: one gene symbol per line,
    trailing '# comment' text stripped, blank/full-line-comment lines
    ignored. Returns a dict mapping UPPERCASE symbol -> original-case
    symbol as first seen (for display), so matching is case-insensitive
    while the report can still show the symbol as the user wrote it."""
    if not Path(path).exists():
        raise FileNotFoundError(
            f"Gene list file not found: {path} - see xenobiotic_resistance_cyp_genes.example.txt "
            "for the format; copy it, rename it, and edit the gene list, then pass it via --gene-list."
        )
    symbols = {}
    with open(path) as fh:
        for line in fh:
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            if line.upper() not in symbols:
                symbols[line.upper()] = line
    if not symbols:
        raise ValueError(f"No gene symbols found in {path} (every line was blank or a comment)")
    return symbols


def restrict_to_gene_list(gene_rows, target_symbols_upper):
    """Filters a species' full per-gene row list down to just the genes
    whose symbol (case-insensitively) is in target_symbols_upper (a set
    of UPPERCASE symbols). Returns (kept_rows, found_upper_set,
    missing_upper_set) - missing is which target genes were NOT found at
    all in this species' known Cyp gene set (naming-mismatch diagnostic)."""
    kept = []
    found_upper = set()
    for g in gene_rows:
        su = g["symbol"].upper()
        if su in target_symbols_upper:
            kept.append(g)
            found_upper.add(su)
    missing_upper = target_symbols_upper - found_upper
    return kept, found_upper, missing_upper


# ==========================================================================
# 2. Report generation

def build_xenobiotic_report(species_summaries, pooled, alpha, target_symbols,
                             missing_by_species, gene_list_path, boot=None):
    lines = []
    lines.append("# TE concentration in literature-curated insecticide-resistance Cyp genes")
    lines.append("")
    lines.append(
        f"This report restricts the Cyp-gene TE comparison to a specific, literature-curated "
        f"list of {len(target_symbols)} gene(s) with established roles in insecticide/xenobiotic "
        f"resistance (loaded from `{gene_list_path}`), rather than every annotated Cyp gene."
    )
    lines.append("")
    lines.append(
        "> **This gene list is a starting point, not an authoritative or exhaustive one** - "
        "review it against your own literature search. See the comments in the gene-list file "
        "for sources/caveats on each included gene."
    )
    lines.append("")
    lines.append("**Target genes:** " + ", ".join(sorted(target_symbols.values(), key=str.lower)))
    lines.append("")

    lines.append("## Per-species summary (target genes only)")
    lines.append("")
    lines.append(
        "| Species | Exposure | Target genes found | Genes w/ TE | % w/ TE | Total TEs | TE/gene |"
    )
    lines.append("|---|---|---:|---:|---:|---:|---:|")
    for s in species_summaries:
        lines.append(
            f"| {s['species']} | {s['exposure']} | {s['n_cyp_genes']} of {len(target_symbols)} | "
            f"{s['n_genes_with_te']} | {s['pct_genes_with_te']:.1f}% | {s['total_te_count']} | "
            f"{s['te_per_gene']:.3f} |"
        )
    lines.append("")

    for sp, missing_upper in missing_by_species.items():
        if missing_upper:
            missing_display = sorted((target_symbols.get(u, u) for u in missing_upper), key=str.lower)
            lines.append(
                f"> Note: {sp} - {len(missing_upper)} target gene(s) not found in this species' "
                f"`gene` features at all: {missing_display}. Likely a naming/ortholog mismatch "
                "(check for paralog-suffixed variants like Cyp12d1-d/Cyp12d1-p) rather than the "
                "gene being genuinely absent from the genome."
            )
    lines.append("")

    if pooled is None:
        lines.append(
            "## Pooled per-gene statistical test\n\n"
            "Not run - fewer than 1 species with usable target-gene data in one of the two "
            "exposure groups (check the per-species table above for 0-gene rows)."
        )
        return "\n".join(lines)

    lines.append("## Pooled per-gene statistical test (target genes only)")
    lines.append("")
    lines.append(
        "Every matched target gene across all species with usable data is pooled into one table "
        f"of {pooled['high_n'] + pooled['low_n']} gene-level observations ({pooled['high_n']} "
        f"high-exposure, {pooled['low_n']} low-exposure):"
    )
    lines.append("")
    t = pooled["table"]
    lines.append("**2x2 table (exposure group x gene has >=1 TE):**")
    lines.append("")
    lines.append("| | Has TE | No TE |")
    lines.append("|---|---:|---:|")
    lines.append(f"| High exposure | {t['high_has_te']} | {t['high_no_te']} |")
    lines.append(f"| Low exposure | {t['low_has_te']} | {t['low_no_te']} |")
    lines.append("")
    lines.append(f"- **Fisher's exact test (two-tailed):** p = {pooled['fisher_p']:.6g}")
    lines.append(f"- Chi-square (Yates-corrected, cross-check): chi2 = {pooled['chi2']:.4f}, p = {pooled['chi2_p']:.6g}")
    lines.append("")
    lines.append(
        f"**Mann-Whitney U test on per-gene TE counts** (high-exposure target genes averaged "
        f"{pooled['high_mean_te']:.3f} TEs/gene vs {pooled['low_mean_te']:.3f} TEs/gene for "
        f"low-exposure): U = {pooled['mwu_U']:.1f}, p = {pooled['mwu_p']:.6g}"
    )
    lines.append(
        f"- **Welch's t-test (unequal variances) on per-gene TE counts:** "
        f"t = {pooled['ttest_t']:.3f}, df = {pooled['ttest_df']:.1f}, p = {pooled['ttest_p']:.6g} "
        "*(a familiar reference point - Mann-Whitney above is the more statistically appropriate "
        "test for this skewed count data, especially with such a small curated gene list)*"
    )
    lines.append("")

    table_saturated = (t["high_no_te"] == 0 and t["low_no_te"] == 0)
    if table_saturated:
        lines.append(
            "> Note: every matched target gene has >=1 associated TE, so the presence/absence "
            "table above has no variance to test - Fisher's/chi-square p=1 is an artifact of "
            "that saturation. The Mann-Whitney U test on TE counts is the test actually carrying "
            "signal here."
        )
        lines.append("")

    fisher_sig = pooled["fisher_p"] < alpha
    mwu_sig = pooled["mwu_p"] < alpha
    if pooled["high_mean_te"] > pooled["low_mean_te"]:
        direction_supports = True
    elif pooled["high_mean_te"] < pooled["low_mean_te"]:
        direction_supports = False
    else:
        direction_supports = None

    lines.append(f"## Verdict (alpha = {alpha})")
    lines.append("")
    if direction_supports is None:
        verdict = (
            "**Inconclusive (tied).** High- and low-exposure target genes show identical mean "
            "TE burden - no directional signal either way."
        )
    elif fisher_sig and mwu_sig and direction_supports:
        verdict = (
            "**Supports the hypothesis, within the curated resistance-gene list.** Both tests "
            "are significant and in the predicted direction: high-exposure species' known "
            "resistance-associated Cyp genes carry more TEs than low-exposure species'."
        )
    elif (fisher_sig or mwu_sig) and direction_supports:
        verdict = (
            "**Partially supports the hypothesis, within the curated resistance-gene list.** "
            "Direction is as predicted and one of the two tests reaches significance, but not "
            "both - suggestive, not conclusive."
            + (" The presence/absence test is uninformative here (saturated table, see note "
               "above), so the Mann-Whitney result carries the real weight." if table_saturated else "")
        )
    elif not direction_supports and (fisher_sig or mwu_sig):
        verdict = (
            "**Contradicts the hypothesis, within the curated resistance-gene list.** At least "
            "one test is significant but in the OPPOSITE direction from predicted."
        )
    else:
        verdict = (
            "**Inconclusive.** Neither test reaches significance within this small, curated gene "
            "set - with only a handful of genes, this test has very little statistical power "
            "either way."
        )
    lines.append(verdict)
    lines.append("")

    if boot is not None:
        lines.append(base.format_bootstrap_section(boot, alpha))

    lines.append("## Caveats")
    lines.append("")
    lines.append(
        "- **Gene list curation:** results are only as good as the target gene list - it's a "
        "literature starting point (see the .example file's citations), not a validated, "
        "exhaustive, or peer-reviewed set for your specific study. Genes with no established "
        "resistance role that you've added, or well-established ones you haven't, will change "
        "these results."
    )
    lines.append(
        "- **Very small sample size:** a curated list is necessarily much smaller than the full "
        "Cyp gene set (often single digits to low tens of genes per species) - both tests have "
        "correspondingly low statistical power, so a non-significant result here is weak "
        "evidence of \"no effect,\" and a significant one deserves scrutiny for whether it's "
        "being driven by just one or two genes/species."
    )
    lines.append(
        "- **Naming mismatches silently reduce the effective list:** a target gene not found in "
        "a species' annotation (see the notes above) is simply excluded from that species' count "
        "- it does not count as \"no TE\" or otherwise bias the comparison, but it does shrink "
        "the effective sample for that species."
    )
    lines.append(
        "- **Same pseudoreplication/correlation-vs-causation caveats as the other comparisons "
        "apply here** (see compare_te_cyp_exposure.py's report) - association shown is not proof "
        "of causation, and pooling genes across species treats them as independent when they "
        "share ancestry."
    )
    lines.append("")

    return "\n".join(lines)


# ==========================================================================
# 3. Main

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Compare TE concentration in a literature-curated list of insecticide-"
                    "resistance Cyp genes across species, using the combined GFF3 outputs from "
                    "build_tfbs_te_gff.py.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--config", help="Path to species config INI (same format/file as compare_te_cyp_exposure.py)")
    parser.add_argument("--gene-list", help="Path to the target gene list text file "
                                             "(see xenobiotic_resistance_cyp_genes.example.txt)")
    parser.add_argument("--output-csv", default="te_cyp_xenobiotic_species_summary.csv",
                         help="Where to write the per-species summary CSV (default: %(default)s)")
    parser.add_argument("--per-gene-csv", default=None,
                         help="Optional: also write a per-gene CSV for the target gene subset")
    parser.add_argument("--output-report", default="te_cyp_xenobiotic_report.md",
                         help="Where to write the markdown report (default: %(default)s)")
    parser.add_argument("--alpha", type=float, default=0.05, help="Significance threshold (default: 0.05)")
    parser.add_argument("--plot", dest="plot", action="store_true", default=True,
                         help="Write a bar chart of TE/gene by species (default: on)")
    parser.add_argument("--no-plot", dest="plot", action="store_false", help="Skip the plot entirely")
    parser.add_argument("--plot-path", default="te_cyp_xenobiotic_by_species.png", help="Plot output path (default: %(default)s)")
    parser.add_argument("--self-test", action="store_true",
                         help="Run the internal statistics cross-validation checks (shared with "
                              "compare_te_cyp_exposure.py) and exit")
    args = parser.parse_args(argv)

    if args.self_test:
        ok = base.validate_stats()
        sys.exit(0 if ok else 1)

    if not args.config:
        parser.error("--config is required (unless using --self-test)")
    if not args.gene_list:
        parser.error("--gene-list is required (unless using --self-test) - see "
                      "xenobiotic_resistance_cyp_genes.example.txt")

    ok = base.validate_stats()
    if not ok:
        print("[error] statistics self-test failed - refusing to report results from untrusted code", file=sys.stderr)
        sys.exit(1)
    print()

    target_symbols = parse_gene_list(args.gene_list)
    target_symbols_upper = set(target_symbols.keys())
    print(f"[info] loaded {len(target_symbols)} target gene(s) from {args.gene_list}")

    species_cfg = base.load_species_config(args.config)

    species_summaries = []
    all_gene_rows = []
    per_gene_csv_rows = []
    missing_by_species = {}

    for species_name, cfg in sorted(species_cfg.items()):
        gff_path = cfg["gff"]
        exposure = cfg["exposure"]
        if not Path(gff_path).exists():
            print(f"[warn] {species_name}: GFF3 file not found, skipping: {gff_path}", file=sys.stderr)
            continue

        gene_name_pattern = cfg.get("gene_name_pattern")
        gene_rows, unmatched_te, n_seen, n_kept = base.parse_species_gff3(gff_path, gene_name_pattern)

        target_gene_rows, found_upper, missing_upper = restrict_to_gene_list(gene_rows, target_symbols_upper)
        missing_by_species[species_name] = missing_upper

        print(f"[info] {species_name}: {len(target_gene_rows)} of {len(target_symbols)} target "
              f"genes found ({len(missing_upper)} missing)")

        summary = base.compute_species_summary(species_name, exposure, target_gene_rows)
        species_summaries.append(summary)
        for g in target_gene_rows:
            all_gene_rows.append((species_name, exposure, g))
            per_gene_csv_rows.append({
                "species": species_name, "exposure": exposure, "symbol": g["symbol"],
                "length_bp": g["length_bp"], "te_count": g["te_count"],
            })

    if not species_summaries:
        print("[error] no species data loaded - check --config paths", file=sys.stderr)
        sys.exit(1)

    n_high = sum(1 for s in species_summaries if s["exposure"] == "high" and s["n_cyp_genes"] > 0)
    n_low = sum(1 for s in species_summaries if s["exposure"] == "low" and s["n_cyp_genes"] > 0)
    pooled = None
    boot = None
    if n_high > 0 and n_low > 0:
        pooled = base.run_pooled_tests(all_gene_rows)
        boot = base.bootstrap_species_cluster_ci(all_gene_rows)
    else:
        print("[warn] need at least one species with >=1 matched target gene in each exposure "
              "group to run the pooled test - skipping statistics", file=sys.stderr)

    csv_fields = ["species", "exposure", "n_cyp_genes", "n_genes_with_te", "pct_genes_with_te",
                  "total_te_count", "total_cyp_bp", "te_per_kb", "te_per_gene"]
    with open(args.output_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=csv_fields)
        writer.writeheader()
        for s in species_summaries:
            writer.writerow({k: s[k] for k in csv_fields})
    print(f"[info] wrote per-species summary CSV to {args.output_csv}")

    if args.per_gene_csv:
        with open(args.per_gene_csv, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=["species", "exposure", "symbol", "length_bp", "te_count"])
            writer.writeheader()
            writer.writerows(per_gene_csv_rows)
        print(f"[info] wrote per-gene CSV to {args.per_gene_csv}")

    report_text = build_xenobiotic_report(species_summaries, pooled, args.alpha, target_symbols,
                                           missing_by_species, args.gene_list, boot)
    with open(args.output_report, "w") as fh:
        fh.write(report_text)
    print(f"[info] wrote report to {args.output_report}")

    if args.plot:
        base.maybe_write_plot(species_summaries, args.plot_path)

    if pooled:
        print()
        print(f"Fisher's exact p = {pooled['fisher_p']:.6g}   Mann-Whitney U p = {pooled['mwu_p']:.6g}")
        print(f"Species-cluster bootstrap: {boot['prob_high_greater']*100:.1f}% of resamples showed "
              f"higher TE/gene in high-exposure species")


if __name__ == "__main__":
    main()
