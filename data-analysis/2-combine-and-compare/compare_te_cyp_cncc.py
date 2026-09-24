#!/usr/bin/env python3
"""
compare_te_cyp_cncc.py

Compares transposable element (TE) counts in Cyp genes ACROSS THE 5
SPECIES, restricted to just the Cyp genes that sit near a predicted
Cap-n-Collar (CncC:Maf-S) antioxidant/xenobiotic response element (ARE) -
the literature-derived motif build_tfbs_te_gff.py adds to its FIMO scan
(see that script's "4b-bis" section) because CncC:Maf-S is the master
regulator of Cyp-mediated xenobiotic/insecticide detoxification in
insects. This narrows the earlier, broader comparison
(compare_te_cyp_exposure.py, which used every Cyp gene) down to the
subset that's plausibly under CncC-driven inducible regulation - the Cyp
genes most directly implicated in the detox response the abstract is
about.

REQUIRES: compare_te_cyp_exposure.py in the same directory - this script
imports its already-validated GFF3 parsing and statistics functions
rather than re-implementing/duplicating them.

INPUT
-----
Same combined GFF3 files (and same --config file) as
compare_te_cyp_exposure.py. In addition to the `gene` and
`mobile_genetic_element`/TE-Description records that script already
reads, this one also needs `TF_binding_site` records - specifically ones
build_tfbs_te_gff.py writes for the CncC:Maf-S ARE motif, identifiable by
either:
    motif_source_id=CncC_Maf_ARE   (the attribute key/value the literature
                                     motif specifically uses - JASPAR-
                                     sourced hits use jaspar_matrix_id
                                     instead, so this is unambiguous)
    Name=<something containing "cnc">   (fallback, in case a future JASPAR
                                          release adds a real Cnc-family
                                          motif)
A gene counts as "close to a cap-n-collar site" if it appears in that
hit's Description="...near <gene1>,<gene2>,..." attribute (same
convention build_tfbs_te_gff.py uses for every TFBS/TE record).

IMPORTANT DATA CAVEAT: a species' GFF3 only has TF_binding_site records
at all if it was run through build_tfbs_te_gff.py's full pipeline WITH
FIMO motif scanning enabled (not --skip-tfbs, and with fimo actually
available). If a species' file has zero TF_binding_site records, this
script reports 0 CncC-associated Cyp genes for it - but that's a DATA GAP
(motif scanning was never run), not a biological finding, and the report
flags it explicitly so it isn't misread as "this species has no CncC
sites near its Cyp genes."

OUTPUT
------
Same shape as compare_te_cyp_exposure.py: a per-species CSV, a markdown
report (with the same pooled Fisher's exact / Mann-Whitney U tests, now
run only on the CncC-associated gene subset), and an optional bar plot.

USAGE
-----
    python compare_te_cyp_cncc.py --config te_cyp_species_config.ini

    python compare_te_cyp_cncc.py --config te_cyp_species_config.ini \\
        --output-csv te_cyp_cncc_summary.csv \\
        --output-report te_cyp_cncc_report.md \\
        --per-gene-csv te_cyp_cncc_per_gene.csv \\
        --plot

    python compare_te_cyp_cncc.py --self-test

Uses the same --config file as compare_te_cyp_exposure.py (see
te_cyp_species_config.example.ini) - no separate config format needed.
"""

import argparse
import csv
import re
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
# 1. CncC (cap-n-collar) TFBS parsing

CNCC_MOTIF_SOURCE_ID = "CncC_Maf_ARE"   # matches build_tfbs_te_gff.py's CNCC_MAF_ARE_MOTIF_ID
_CNCC_NAME_FALLBACK_RE = re.compile(r"\bcnc\b", re.IGNORECASE)
_NEAR_GENES_RE = re.compile(r"near ([^\s;]+)")


def parse_cncc_associated_genes(path):
    """Scans a species' combined GFF3 for TF_binding_site records that are
    CncC:Maf-S ARE hits (identified by motif_source_id=CncC_Maf_ARE, or a
    Name containing "cnc" as a fallback), and returns the set of gene
    symbols named in their Description="...near <genes>" attribute -
    i.e. the Cyp genes considered "close to a cap-n-collar site".

    Returns (cncc_genes: set[str], n_tfbs_total: int, n_cncc_hits: int).
    n_tfbs_total lets the caller detect the "this species has zero
    TF_binding_site records at all" data-gap case and flag it distinctly
    from "this species genuinely has no CncC-proximal Cyp genes".
    """
    cncc_genes = set()
    n_tfbs_total = 0
    n_cncc_hits = 0

    for seqid, source, ftype, start, end, score, strand, phase, attrs in base._gff3_feature_iter(path):
        if ftype != "TF_binding_site":
            continue
        n_tfbs_total += 1

        motif_source_id = attrs.get("motif_source_id", "")
        name = attrs.get("Name", "")
        is_cncc = (motif_source_id == CNCC_MOTIF_SOURCE_ID) or bool(_CNCC_NAME_FALLBACK_RE.search(name))
        if not is_cncc:
            continue
        n_cncc_hits += 1

        desc = attrs.get("Description", "")
        m = _NEAR_GENES_RE.search(desc)
        if not m:
            continue
        genes_str = m.group(1)
        if not genes_str or genes_str == "NA":
            continue
        for g in genes_str.split(","):
            g = g.strip()
            if g:
                cncc_genes.add(g)

    return cncc_genes, n_tfbs_total, n_cncc_hits


def restrict_to_cncc_genes(gene_rows, cncc_genes):
    """Filters a species' full per-gene row list (from
    compare_te_cyp_exposure.parse_species_gff3) down to just the genes
    that are also in cncc_genes, and reports any cncc_genes entries that
    didn't match a known gene (naming mismatches, same idea as the
    TE-to-gene "unmatched" caveat in the base script)."""
    known_symbols = {g["symbol"] for g in gene_rows}
    kept = [g for g in gene_rows if g["symbol"] in cncc_genes]
    unmatched = sorted(cncc_genes - known_symbols)
    return kept, unmatched


# ==========================================================================
# 2. Report generation (mirrors build_report's structure/tests, but framed
#    around the CncC-restricted subset and the TFBS-coverage caveat)

def build_cncc_report(species_summaries, pooled, alpha, unmatched_cncc_by_species,
                       unmatched_te_by_species, coverage_by_species, boot=None):
    lines = []
    lines.append("# TE counts in CncC (cap-n-collar)-proximal Cyp genes, across species")
    lines.append("")
    lines.append(
        "This report narrows the broader Cyp-gene TE comparison down to just the Cyp genes "
        "that sit near a predicted CncC:Maf-S antioxidant/xenobiotic response element (ARE) - "
        "the master regulator of Cyp-mediated insecticide/xenobiotic detoxification in insects "
        "(Veraksa 2000; Misra et al. 2011). These are the Cyp genes most plausibly under "
        "CncC-driven inducible regulation, i.e. the ones most directly relevant to the "
        "TE-affects-detox-gene-expression hypothesis, rather than every annotated Cyp gene."
    )
    lines.append("")

    # --- coverage caveat: species with zero TF_binding_site records at all ---
    zero_coverage = sorted(sp for sp, n in coverage_by_species.items() if n == 0)
    if zero_coverage:
        verb = "has" if len(zero_coverage) == 1 else "have"
        plural = "" if len(zero_coverage) == 1 else "s"
        lines.append(
            f"> **Data gap, not biology:** {', '.join(zero_coverage)} {verb} zero "
            f"`TF_binding_site` records in {'its' if len(zero_coverage) == 1 else 'their'} GFF3 "
            f"at all - meaning that species' file{plural} {'was' if len(zero_coverage) == 1 else 'were'} "
            "never run through build_tfbs_te_gff.py's FIMO motif-scanning step (not that "
            f"{'it lacks' if len(zero_coverage) == 1 else 'they lack'} CncC sites). "
            f"{'Its' if len(zero_coverage) == 1 else 'Their'} 0 CncC-associated Cyp genes below "
            "should be read as **missing data**, and excluded from any real interpretation, until "
            f"{'that species is' if len(zero_coverage) == 1 else 'those species are'} re-run with "
            "TFBS scanning enabled."
        )
        lines.append("")

    lines.append("## Per-species summary (CncC-proximal Cyp genes only)")
    lines.append("")
    lines.append(
        "| Species | Exposure | TF_binding_site records | CncC-proximal Cyp genes | Genes w/ TE | "
        "% w/ TE | Total TEs | TE/gene |"
    )
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|")
    for s in species_summaries:
        n_tfbs = coverage_by_species.get(s["species"], 0)
        flag = " *(no TFBS data)*" if n_tfbs == 0 else ""
        lines.append(
            f"| {s['species']}{flag} | {s['exposure']} | {n_tfbs} | {s['n_cyp_genes']} | "
            f"{s['n_genes_with_te']} | {s['pct_genes_with_te']:.1f}% | {s['total_te_count']} | "
            f"{s['te_per_gene']:.3f} |"
        )
    lines.append("")

    for sp, unmatched in unmatched_cncc_by_species.items():
        if unmatched:
            lines.append(
                f"> Note: {sp} - {len(unmatched)} CncC-hit gene reference(s) did not match any "
                f"known Cyp gene in that species' `gene` features ({unmatched[:5]}"
                f"{' ...' if len(unmatched) > 5 else ''}); likely a gene-naming mismatch "
                "(e.g. paralog suffixes) rather than a real gap."
            )
    for sp, unmatched in unmatched_te_by_species.items():
        if unmatched:
            lines.append(
                f"> Note: {sp} - {len(unmatched)} TE record(s) referenced a gene symbol not found "
                f"among that species' `gene` features ({unmatched[:5]}"
                f"{' ...' if len(unmatched) > 5 else ''}); excluded from TE counts above."
            )
    lines.append("")

    if pooled is None:
        lines.append(
            "## Pooled per-gene statistical test\n\n"
            "Not run - fewer than 1 species with usable data in one of the two exposure groups "
            "once restricted to CncC-proximal genes (often because of the TF_binding_site data "
            "gap noted above)."
        )
        return "\n".join(lines)

    lines.append("## Pooled per-gene statistical test (CncC-proximal genes only)")
    lines.append("")
    lines.append(
        "Every CncC-proximal Cyp gene across the species with usable TFBS data is pooled into "
        f"one table of {pooled['high_n'] + pooled['low_n']} gene-level observations "
        f"({pooled['high_n']} high-exposure, {pooled['low_n']} low-exposure):"
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
        f"**Mann-Whitney U test on per-gene TE counts** (high-exposure CncC-proximal genes "
        f"averaged {pooled['high_mean_te']:.3f} TEs/gene vs {pooled['low_mean_te']:.3f} TEs/gene "
        f"for low-exposure): U = {pooled['mwu_U']:.1f}, p = {pooled['mwu_p']:.6g}"
    )
    lines.append(
        f"- **Welch's t-test (unequal variances) on per-gene TE counts:** "
        f"t = {pooled['ttest_t']:.3f}, df = {pooled['ttest_df']:.1f}, p = {pooled['ttest_p']:.6g} "
        "*(a familiar reference point - Mann-Whitney above is the more statistically appropriate "
        "test for this skewed count data)*"
    )
    lines.append("")

    table_saturated = (t["high_no_te"] == 0 and t["low_no_te"] == 0)
    if table_saturated:
        lines.append(
            "> Note: every CncC-proximal Cyp gene has >=1 associated TE, so the presence/absence "
            "table above has no variance to test - Fisher's/chi-square p=1 is an artifact of that "
            "saturation, not evidence of no effect. The Mann-Whitney U test on TE counts is the "
            "test actually carrying signal here."
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
            "**Inconclusive (tied).** High- and low-exposure CncC-proximal genes show identical "
            "mean TE burden - no directional signal either way."
        )
    elif fisher_sig and mwu_sig and direction_supports:
        verdict = (
            "**Supports the hypothesis, within the CncC-proximal subset.** Both tests are "
            "significant and in the predicted direction: high-exposure species' CncC-regulated "
            "Cyp genes carry more TEs than low-exposure species'."
        )
    elif (fisher_sig or mwu_sig) and direction_supports:
        verdict = (
            "**Partially supports the hypothesis, within the CncC-proximal subset.** Direction is "
            "as predicted and one of the two tests reaches significance, but not both - "
            "suggestive, not conclusive."
            + (" The presence/absence test is uninformative here (saturated table, see note "
               "above), so the Mann-Whitney result carries the real weight." if table_saturated else "")
        )
    elif not direction_supports and (fisher_sig or mwu_sig):
        verdict = (
            "**Contradicts the hypothesis, within the CncC-proximal subset.** At least one test "
            "is significant but in the OPPOSITE direction from predicted."
        )
    else:
        verdict = (
            "**Inconclusive.** Neither test reaches significance within this (typically much "
            "smaller) CncC-proximal gene subset."
        )
    lines.append(verdict)
    lines.append("")

    if boot is not None:
        lines.append(base.format_bootstrap_section(boot, alpha))

    lines.append("## Caveats")
    lines.append("")
    lines.append(
        "- **Same pseudoreplication/correlation-vs-causation caveats as the full-Cyp-gene "
        "comparison apply here** (see compare_te_cyp_exposure.py's report) - pooling genes "
        "across species treats them as independent when they share ancestry, and any "
        "association shown is not proof of causation."
    )
    lines.append(
        "- **Sample size:** the CncC-proximal subset is, by definition, smaller than the full "
        "Cyp gene set - with fewer genes, both tests have less power, so a non-significant "
        "result here is weaker evidence of \"no effect\" than the same result would be on the "
        "full gene set."
    )
    lines.append(
        "- **CncC-association source:** as with TE-to-gene assignment, which genes count as "
        "\"near\" a CncC site is inherited as-is from build_tfbs_te_gff.py's FIMO scan window "
        "and promoter-region definition (--upstream/--downstream/--te-flank) - this script does "
        "not re-derive or second-guess that call."
    )
    lines.append(
        "- **Motif is literature-derived, not empirical:** the CncC:Maf-S ARE motif is a "
        "consensus PFM built from published sequence data, not from ChIP-seq/SELEX like the "
        "JASPAR motifs - treat CncC \"hits\" as candidate sites worth follow-up, not "
        "experimentally confirmed binding events."
    )
    lines.append("")

    return "\n".join(lines)


# ==========================================================================
# 3. Main

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Compare TE counts in CncC (cap-n-collar)-proximal Cyp genes across species, "
                    "using the combined GFF3 outputs from build_tfbs_te_gff.py.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--config", help="Path to species config INI (same format/file as compare_te_cyp_exposure.py)")
    parser.add_argument("--output-csv", default="te_cyp_cncc_species_summary.csv",
                         help="Where to write the per-species summary CSV (default: %(default)s)")
    parser.add_argument("--per-gene-csv", default=None,
                         help="Optional: also write a per-gene CSV for the CncC-proximal subset")
    parser.add_argument("--output-report", default="te_cyp_cncc_report.md",
                         help="Where to write the markdown report (default: %(default)s)")
    parser.add_argument("--alpha", type=float, default=0.05, help="Significance threshold (default: 0.05)")
    parser.add_argument("--plot", dest="plot", action="store_true", default=True,
                         help="Write a bar chart of TE/gene by species (default: on)")
    parser.add_argument("--no-plot", dest="plot", action="store_false", help="Skip the plot entirely")
    parser.add_argument("--plot-path", default="te_cyp_cncc_by_species.png", help="Plot output path (default: %(default)s)")
    parser.add_argument("--self-test", action="store_true",
                         help="Run the internal statistics cross-validation checks (shared with "
                              "compare_te_cyp_exposure.py) and exit")
    args = parser.parse_args(argv)

    if args.self_test:
        ok = base.validate_stats()
        sys.exit(0 if ok else 1)

    if not args.config:
        parser.error("--config is required (unless using --self-test)")

    ok = base.validate_stats()
    if not ok:
        print("[error] statistics self-test failed - refusing to report results from untrusted code", file=sys.stderr)
        sys.exit(1)
    print()

    species_cfg = base.load_species_config(args.config)

    species_summaries = []
    all_gene_rows = []          # (species, exposure, gene_row_dict) - CncC-proximal subset only
    per_gene_csv_rows = []
    unmatched_cncc_by_species = {}
    unmatched_te_by_species = {}
    coverage_by_species = {}    # species -> n_tfbs_total (for the data-gap flag)

    for species_name, cfg in sorted(species_cfg.items()):
        gff_path = cfg["gff"]
        exposure = cfg["exposure"]
        if not Path(gff_path).exists():
            print(f"[warn] {species_name}: GFF3 file not found, skipping: {gff_path}", file=sys.stderr)
            continue

        gene_name_pattern = cfg.get("gene_name_pattern")
        gene_rows, unmatched_te, n_seen, n_kept = base.parse_species_gff3(gff_path, gene_name_pattern)
        unmatched_te_by_species[species_name] = unmatched_te

        cncc_genes, n_tfbs_total, n_cncc_hits = parse_cncc_associated_genes(gff_path)
        coverage_by_species[species_name] = n_tfbs_total

        cncc_gene_rows, unmatched_cncc = restrict_to_cncc_genes(gene_rows, cncc_genes)
        unmatched_cncc_by_species[species_name] = unmatched_cncc

        if n_tfbs_total == 0:
            print(f"[warn] {species_name}: 0 TF_binding_site records found in {gff_path} - this "
                  f"species was likely never run through FIMO motif scanning; CncC-proximal gene "
                  f"count below is a DATA GAP, not a real zero", file=sys.stderr)
        else:
            print(f"[info] {species_name}: {n_tfbs_total} TF_binding_site records, "
                  f"{n_cncc_hits} CncC:Maf-S ARE hits, {len(cncc_gene_rows)} of {len(gene_rows)} "
                  f"Cyp genes are CncC-proximal")

        summary = base.compute_species_summary(species_name, exposure, cncc_gene_rows)
        species_summaries.append(summary)
        for g in cncc_gene_rows:
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
        print("[warn] need at least one species with >=1 CncC-proximal gene in each exposure "
              "group to run the pooled test - skipping statistics (see report for why)", file=sys.stderr)

    # --- per-species CSV ---
    csv_fields = ["species", "exposure", "n_tfbs_total", "n_cyp_genes", "n_genes_with_te",
                  "pct_genes_with_te", "total_te_count", "total_cyp_bp", "te_per_kb", "te_per_gene"]
    with open(args.output_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=csv_fields)
        writer.writeheader()
        for s in species_summaries:
            row = {k: s[k] for k in csv_fields if k in s}
            row["n_tfbs_total"] = coverage_by_species.get(s["species"], 0)
            writer.writerow(row)
    print(f"[info] wrote per-species summary CSV to {args.output_csv}")

    if args.per_gene_csv:
        with open(args.per_gene_csv, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=["species", "exposure", "symbol", "length_bp", "te_count"])
            writer.writeheader()
            writer.writerows(per_gene_csv_rows)
        print(f"[info] wrote per-gene CSV to {args.per_gene_csv}")

    # --- report ---
    report_text = build_cncc_report(species_summaries, pooled, args.alpha,
                                     unmatched_cncc_by_species, unmatched_te_by_species,
                                     coverage_by_species, boot)
    with open(args.output_report, "w") as fh:
        fh.write(report_text)
    print(f"[info] wrote report to {args.output_report}")

    # --- plot ---
    if args.plot:
        base.maybe_write_plot(species_summaries, args.plot_path)

    if pooled:
        print()
        print(f"Fisher's exact p = {pooled['fisher_p']:.6g}   Mann-Whitney U p = {pooled['mwu_p']:.6g}")
        print(f"Species-cluster bootstrap: {boot['prob_high_greater']*100:.1f}% of resamples showed "
              f"higher TE/gene in high-exposure species")


if __name__ == "__main__":
    main()
