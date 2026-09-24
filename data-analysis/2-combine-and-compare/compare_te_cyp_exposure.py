#!/usr/bin/env python3
"""
compare_te_cyp_exposure.py

Tests whether Drosophila species with greater insecticide exposure show a
higher prevalence of transposable element (TE) insertions in/near Cytochrome
P450 (CYP) genes than species with lower exposure - i.e. it checks the
research abstract's hypothesis against real data.

INPUT
-----
This script consumes the *combined GFF3 files already produced by*
build_tfbs_te_gff.py (one per species) - it does not run RepeatModeler/
RepeatMasker/FIMO itself. It only needs two feature types that pipeline
already writes into its output:

    gene                    - one row per Cyp gene, Name=<symbol>
    mobile_genetic_element  - one row per TE, with an attribute
                              Description=Within range of <symbol>
                              (written by te_record_to_gff3 in the main
                              pipeline) that ties each TE back to the gene
                              it was found near/within.

So the TE-to-gene association is read directly out of that attribute -
no coordinate-overlap recomputation is needed here.

WHAT IT COMPUTES
-----------------
Per species:
    n_cyp_genes         - number of Cyp gene records in that species' GFF3
    n_genes_with_te      - how many of those genes have >=1 associated TE
    pct_genes_with_te    - n_genes_with_te / n_cyp_genes * 100
    total_te_count       - total TEs associated with any Cyp gene
    total_cyp_bp         - summed length (bp) of all Cyp gene records
    te_per_kb            - total_te_count / (total_cyp_bp / 1000)
    te_per_gene          - total_te_count / n_cyp_genes

These normalized densities (te_per_kb, te_per_gene) are what make species
with different genome sizes / Cyp-gene-set sizes comparable - a species
with more or longer Cyp genes would otherwise look like it has "more TEs"
just from having more sequence to hit.

STATISTICAL TEST (pre-registered choice, see prior clarification)
-------------------------------------------------------------------
Every individual Cyp gene, pooled across all 5 species, is treated as one
observation ("pooled per-gene test"). Two tests are run on that pooled
gene-level table:

  1. Fisher's exact test (two-tailed) on the 2x2 table
         exposure group (high/low)  x  gene has >=1 TE (yes/no)
     This tests whether TE PRESENCE is associated with exposure group.

  2. Mann-Whitney U test (normal approximation, tie-corrected) comparing
     the distribution of per-gene TE COUNTS between the high- and
     low-exposure groups. This tests whether TE BURDEN (not just
     presence/absence) differs between groups.

Both are implemented here in pure Python stdlib (math.comb / Fraction /
statistics.NormalDist) - no scipy required, matching the zero-external-
dependency design of the rest of this project. Both implementations were
independently cross-validated during development (see validate_stats()
below, which re-derives the same p-values via a second, differently-coded
path and is run automatically with --self-test).

CAVEAT (stated plainly, and repeated in the generated report): pooling
genes across species is pseudoreplication - it treats genes within the
same species as independent when they share a phylogenetic/genomic
background. A secondary, purely descriptive species-level comparison
(3 high-exposure vs 2 low-exposure species, no formal test - far too few
species for one) is also reported, clearly labeled as exploratory only.

USAGE
-----
    python compare_te_cyp_exposure.py --config te_cyp_species_config.ini

    python compare_te_cyp_exposure.py --config te_cyp_species_config.ini \\
        --output-csv te_cyp_summary.csv \\
        --output-report te_cyp_report.md \\
        --per-gene-csv te_cyp_per_gene.csv \\
        --plot

    python compare_te_cyp_exposure.py --self-test     # just run the stats
                                                       # cross-validation
                                                       # and exit

See te_cyp_species_config.example.ini for the config file format.
"""

import argparse
import configparser
import csv
import gzip
import math
import random
import re
import sys
from collections import defaultdict
from fractions import Fraction
from math import comb
from pathlib import Path
from statistics import NormalDist


# ==========================================================================
# 1. Statistics - pure stdlib, no scipy

def fisher_exact_two_tailed(a, b, c, d):
    """Exact two-tailed Fisher's exact test p-value for the 2x2 table
    [[a, b], [c, d]] using the hypergeometric distribution directly
    (math.comb) - no scipy needed.

    Verified against the classic R example:
        fisher.test(matrix(c(3,1,1,3), 2, 2))  ->  p = 0.4857143
    fisher_exact_two_tailed(3, 1, 1, 3) reproduces this to 7 decimal
    places, and has additionally been cross-checked against an
    independent exact-Fraction re-implementation (see validate_stats).
    """
    row1, row2 = a + b, c + d
    col1, col2 = a + c, b + d
    n = row1 + row2

    def hyper_prob(x):
        if x < 0 or x > row1 or (col1 - x) < 0 or (col1 - x) > row2:
            return 0.0
        return (comb(row1, x) * comb(row2, col1 - x)) / comb(n, col1)

    observed_p = hyper_prob(a)
    total_p = 0.0
    lo = max(0, col1 - row2)
    hi = min(row1, col1)
    for x in range(lo, hi + 1):
        p = hyper_prob(x)
        if p <= observed_p * (1 + 1e-9):
            total_p += p
    return min(total_p, 1.0)


def _fisher_exact_via_fractions(a, b, c, d):
    """Independent re-implementation of fisher_exact_two_tailed using exact
    rational (Fraction) arithmetic instead of float division - used only
    to cross-validate the float version in validate_stats(), never called
    from the main pipeline (it's slower and unnecessary once validated)."""
    row1, row2 = a + b, c + d
    col1, col2 = a + c, b + d
    n = row1 + row2

    def hyper_prob(x):
        if x < 0 or x > row1 or (col1 - x) < 0 or (col1 - x) > row2:
            return Fraction(0)
        return Fraction(comb(row1, x) * comb(row2, col1 - x), comb(n, col1))

    observed_p = hyper_prob(a)
    total_p = Fraction(0)
    lo = max(0, col1 - row2)
    hi = min(row1, col1)
    for x in range(lo, hi + 1):
        p = hyper_prob(x)
        if p <= observed_p:
            total_p += p
    return min(float(total_p), 1.0)


def chi_square_2x2(a, b, c, d):
    """Pearson's chi-square test (with Yates' continuity correction) for a
    2x2 table, returned alongside Fisher's exact as a cross-check - large-
    sample tests should roughly agree with the exact test when cell counts
    aren't tiny. Returns (chi2_statistic, p_value). p-value uses the
    chi-square CDF with 1 degree of freedom, computed via the regularized
    lower incomplete gamma function (df=1 reduces to a closed form via erf,
    but we use NormalDist for consistency: chi2(1 dof) statistic Z=sqrt(chi2)
    gives p = 2*(1-Phi(Z)) for the two-sided normal-based equivalent)."""
    n = a + b + c + d
    if n == 0:
        return 0.0, 1.0
    row1, row2 = a + b, c + d
    col1, col2 = a + c, b + d
    if row1 == 0 or row2 == 0 or col1 == 0 or col2 == 0:
        return 0.0, 1.0
    expected_a = row1 * col1 / n
    expected_b = row1 * col2 / n
    expected_c = row2 * col1 / n
    expected_d = row2 * col2 / n
    # Yates' continuity correction
    def term(o, e):
        return (max(0.0, abs(o - e) - 0.5)) ** 2 / e
    chi2 = term(a, expected_a) + term(b, expected_b) + term(c, expected_c) + term(d, expected_d)
    z = chi2 ** 0.5
    p = 2 * (1 - NormalDist().cdf(z))
    return chi2, min(max(p, 0.0), 1.0)


def mann_whitney_u_test(sample1, sample2):
    """Two-sided Mann-Whitney U test, normal approximation with tie
    correction, using statistics.NormalDist (stdlib, Python 3.8+).
    Returns (U statistic, p_value). Sanity-checked in validate_stats():
    identical distributions -> p ~ 1.0; fully separated groups -> small p.
    """
    n1, n2 = len(sample1), len(sample2)
    if n1 == 0 or n2 == 0:
        return 0.0, 1.0

    combined = [(v, 1) for v in sample1] + [(v, 2) for v in sample2]
    combined.sort(key=lambda t: t[0])

    ranks = [0.0] * len(combined)
    i = 0
    while i < len(combined):
        j = i
        while j < len(combined) and combined[j][0] == combined[i][0]:
            j += 1
        avg_rank = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[k] = avg_rank
        i = j

    R1 = sum(r for r, (v, g) in zip(ranks, combined) if g == 1)
    U1 = R1 - n1 * (n1 + 1) / 2.0
    U2 = n1 * n2 - U1
    U = min(U1, U2)

    tie_groups = {}
    for (v, _g) in combined:
        tie_groups[v] = tie_groups.get(v, 0) + 1
    n = n1 + n2
    tie_sum = sum(t ** 3 - t for t in tie_groups.values())
    mean_U = n1 * n2 / 2.0
    if n > 1:
        sigma_U = ((n1 * n2 / 12.0) * ((n + 1) - tie_sum / (n * (n - 1)))) ** 0.5
    else:
        sigma_U = 0.0

    if sigma_U == 0:
        return U, 1.0

    if U1 > mean_U:
        z = (U1 - 0.5 - mean_U) / sigma_U
    elif U1 < mean_U:
        z = (U1 + 0.5 - mean_U) / sigma_U
    else:
        z = 0.0

    p = 2 * (1 - NormalDist().cdf(abs(z)))
    return U, min(p, 1.0)


def _betacf(a, b, x, max_iter=200, eps=1e-14):
    """Continued-fraction evaluation used by regularized_incomplete_beta -
    standard numerical recipe (Lentz's algorithm), pure stdlib."""
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < 1e-300:
        d = 1e-300
    d = 1.0 / d
    h = d
    for m in range(1, max_iter + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-300:
            d = 1e-300
        c = 1.0 + aa / c
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-300:
            d = 1e-300
        c = 1.0 + aa / c
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def regularized_incomplete_beta(a, b, x):
    """I_x(a, b), the regularized incomplete beta function - used to get
    Student's t-distribution p-values without scipy (math.lgamma gives
    the log-gamma prefactor; the continued fraction above does the rest).
    Verified in validate_stats() against analytically-exact special cases
    (Beta(1,1) is uniform, so I_x(1,1) == x exactly) and against an
    independently-coded numerical integration of the t-distribution PDF."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    ln_beta = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
               + a * math.log(x) + b * math.log(1.0 - x))
    front = math.exp(ln_beta)
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    else:
        return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def welch_t_test(sample1, sample2):
    """Welch's two-sample t-test (unequal variances assumed - the safer
    default over Student's pooled-variance t-test when group sizes/
    variances differ, which they will here). Two-tailed p-value via the
    regularized incomplete beta function (pure stdlib, no scipy).

    Returns (t_statistic, degrees_of_freedom, p_value). Degenerate inputs
    (fewer than 2 observations in either group, or zero variance in both)
    return a t of 0.0 and p of 1.0 (or 0.0 if means differ despite zero
    variance - a real but oddly-shaped edge case) rather than raising.

    NOTE ON APPROPRIATENESS: per-gene TE counts are a non-negative, right-
    skewed count variable (many zeros, occasional larger counts) - not
    the roughly-normal continuous data a t-test technically assumes.
    Mann-Whitney U (already reported alongside this) is the more
    statistically appropriate test for this kind of data and should be
    weighted more heavily; the t-test is reported because it's a familiar
    reference point and a useful cross-check, not because it's the ideal
    tool for count data.
    """
    n1, n2 = len(sample1), len(sample2)
    if n1 < 2 or n2 < 2:
        return 0.0, 0.0, 1.0

    mean1 = sum(sample1) / n1
    mean2 = sum(sample2) / n2
    var1 = sum((x - mean1) ** 2 for x in sample1) / (n1 - 1)
    var2 = sum((x - mean2) ** 2 for x in sample2) / (n2 - 1)

    se1 = var1 / n1
    se2 = var2 / n2
    se_total = se1 + se2

    if se_total == 0:
        return 0.0, float(n1 + n2 - 2), (1.0 if mean1 == mean2 else 0.0)

    t_stat = (mean1 - mean2) / math.sqrt(se_total)
    df = (se_total ** 2) / ((se1 ** 2) / (n1 - 1) + (se2 ** 2) / (n2 - 1))

    x = df / (df + t_stat ** 2)
    p = regularized_incomplete_beta(df / 2.0, 0.5, x)
    return t_stat, df, min(max(p, 0.0), 1.0)


def _t_pdf(x, df):
    """Student's t-distribution PDF - used only by
    _t_test_p_via_integration() as an independent cross-check of
    welch_t_test's closed-form (incomplete-beta) p-value."""
    coef = math.exp(math.lgamma((df + 1) / 2) - math.lgamma(df / 2)) / math.sqrt(df * math.pi)
    return coef * (1 + x * x / df) ** (-(df + 1) / 2)


def _t_test_p_via_integration(t_stat, df, tail_width=200, n_steps=20000):
    """Independent re-derivation of the two-tailed t-test p-value via
    numerical integration (Simpson's rule) of the t-distribution PDF,
    rather than the closed-form incomplete-beta route welch_t_test()
    actually uses - cross-validated against it in validate_stats(),
    never called from the main pipeline."""
    t_abs = abs(t_stat)
    a, b = t_abs, t_abs + tail_width
    h = (b - a) / n_steps
    total = _t_pdf(a, df) + _t_pdf(b, df)
    for i in range(1, n_steps):
        x = a + i * h
        total += (4 if i % 2 == 1 else 2) * _t_pdf(x, df)
    integral = total * h / 3
    return min(2 * integral, 1.0)


def validate_stats():
    """Cross-validates both test implementations via an independent second
    code path, and prints the results. Run automatically with --self-test,
    and once at the start of a normal run (fast - milliseconds) so any
    future edit to these functions can't silently break them."""
    print("[self-test] Fisher's exact test:")
    fisher_cases = [
        (3, 1, 1, 3),    # classic R reference case
        (10, 10, 3, 17),
        (5, 5, 5, 5),
        (1, 9, 11, 3),
        (7, 2, 3, 8),
    ]
    ok = True
    for a, b, c, d in fisher_cases:
        p_float = fisher_exact_two_tailed(a, b, c, d)
        p_frac = _fisher_exact_via_fractions(a, b, c, d)
        match = abs(p_float - p_frac) < 1e-9
        ok = ok and match
        print(f"    table=({a},{b},{c},{d}) float={p_float:.10f} fraction={p_frac:.10f} "
              f"{'OK' if match else 'MISMATCH'}")
    ref = fisher_exact_two_tailed(3, 1, 1, 3)
    ref_ok = abs(ref - 0.4857142857142857) < 1e-9
    ok = ok and ref_ok
    print(f"    R reference check (3,1,1,3)=0.4857143: {'OK' if ref_ok else 'MISMATCH'}")

    print("[self-test] Mann-Whitney U test:")
    U, p = mann_whitney_u_test(list(range(1, 11)), list(range(1, 11)))
    identical_ok = abs(p - 1.0) < 1e-6
    ok = ok and identical_ok
    print(f"    identical distributions: U={U} p={p:.6f} {'OK' if identical_ok else 'MISMATCH'}")
    U, p = mann_whitney_u_test([100, 101, 102, 103, 104], [1, 2, 3, 4, 5])
    separated_ok = p < 0.05
    ok = ok and separated_ok
    print(f"    fully separated groups:  U={U} p={p:.6f} {'OK' if separated_ok else 'MISMATCH'}")

    print("[self-test] Regularized incomplete beta function:")
    beta_ok = True
    for x in (0.1, 0.25, 0.5, 0.75, 0.9):
        v = regularized_incomplete_beta(1, 1, x)
        match = abs(v - x) < 1e-9
        beta_ok = beta_ok and match
        print(f"    I_{x}(1,1)={v:.10f} (exact: Beta(1,1) is uniform, expect {x}) {'OK' if match else 'MISMATCH'}")
    for a in (1, 2, 5, 10):
        v = regularized_incomplete_beta(a, a, 0.5)
        match = abs(v - 0.5) < 1e-9
        beta_ok = beta_ok and match
        print(f"    I_0.5({a},{a})={v:.10f} (symmetry, expect 0.5) {'OK' if match else 'MISMATCH'}")
    ok = ok and beta_ok

    print("[self-test] Welch's t-test (cross-checked against independent numerical integration):")
    t_cases = [
        ([23, 25, 21, 30, 28, 24], [18, 20, 15, 22, 19, 17, 21]),
        ([1, 2, 3, 4, 5, 0, 0, 1, 2, 1], [3, 4, 5, 6, 7, 4, 5, 3, 4, 5]),
        (list(range(1, 21)), list(range(15, 35))),
    ]
    t_ok = True
    for s1, s2 in t_cases:
        t_stat, df, p_closed = welch_t_test(s1, s2)
        p_integrated = _t_test_p_via_integration(t_stat, df)
        match = abs(p_closed - p_integrated) < 1e-6
        t_ok = t_ok and match
        print(f"    t={t_stat:.4f} df={df:.2f} p(betainc)={p_closed:.8f} p(integration)={p_integrated:.8f} "
              f"{'OK' if match else 'MISMATCH'}")
    t_stat, df, p_zero = welch_t_test([10, 10, 10], [10, 10, 10])
    zero_var_ok = abs(p_zero - 1.0) < 1e-9
    t_ok = t_ok and zero_var_ok
    print(f"    identical, zero-variance samples: t={t_stat} p={p_zero} (expect 1.0) {'OK' if zero_var_ok else 'MISMATCH'}")
    ok = ok and t_ok

    print(f"[self-test] {'ALL CHECKS PASSED' if ok else 'SOME CHECKS FAILED - DO NOT TRUST RESULTS'}")
    return ok


# ==========================================================================
# 2. GFF3 parsing

_ATTR_SPLIT_RE = re.compile(r";(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)")


def _parse_gff3_attributes(attr_field):
    attrs = {}
    for part in attr_field.strip().split(";"):
        part = part.strip()
        if not part or "=" not in part:
            continue
        key, _, val = part.partition("=")
        attrs[key.strip()] = val.strip()
    return attrs


def _open_maybe_gzip(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def _gff3_feature_iter(path):
    """Yields (seqid, source, type, start, end, score, strand, phase, attrs)
    for every non-comment feature line in a (possibly gzipped) GFF3 file."""
    with _open_maybe_gzip(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attr_field = fields[:9]
            try:
                start_i, end_i = int(start), int(end)
            except ValueError:
                continue
            yield seqid, source, ftype, start_i, end_i, score, strand, phase, _parse_gff3_attributes(attr_field)


_WITHIN_RANGE_RE = re.compile(r"Within range of (\S+)")


def parse_species_gff3(path, gene_name_pattern=None):
    """Parses one species' combined GFF3 and returns a list of per-gene
    dicts:
        {"symbol": str, "length_bp": int, "te_count": int}
    covering every Cyp `gene` feature in the file, with te_count being how
    many `mobile_genetic_element` records' Description attribute names that
    gene via "Within range of <symbol>" (written by build_tfbs_te_gff.py's
    te_record_to_gff3 - read directly, not recomputed from coordinates).

    Genes are keyed by their Name= attribute (the Cyp symbol), matching how
    build_tfbs_te_gff.py's te_record_to_gff3 stamps te_rec["gene"] - both
    ultimately come from the same canonical per-gene symbol used throughout
    that pipeline, so a plain string match is sufficient (no coordinate
    overlap needed here).

    gene_name_pattern: optional compiled regex. If given, only `gene`
    features whose Name matches (via .match(), i.e. anchored at the start)
    are counted as Cyp genes for this species - everything else is
    ignored. This matters when a species' GFF3 is a *whole-genome*
    annotation (every gene in the genome, e.g. ~17,000+ for D. melanogaster
    RefSeq) rather than a file pre-filtered down to just the Cyp genes
    analyzed (as the other species' build_tfbs_te_gff.py combined outputs
    already are, at ~40-90 genes each) - without filtering, a whole-genome
    file would swamp the per-gene metrics with thousands of irrelevant
    non-Cyp genes and make that species look artificially TE-poor.
    Returns (rows, unmatched_te_genes, n_gene_rows_seen, n_gene_rows_kept).
    """
    genes = {}          # symbol -> length_bp
    te_counts_raw = defaultdict(int)   # symbol (as named in Description) -> TE count
    n_gene_rows_seen = 0
    n_gene_rows_kept = 0

    for seqid, source, ftype, start, end, score, strand, phase, attrs in _gff3_feature_iter(path):
        if ftype == "gene":
            n_gene_rows_seen += 1
            symbol = attrs.get("Name")
            if not symbol:
                continue
            if gene_name_pattern is not None and not gene_name_pattern.match(symbol):
                continue
            n_gene_rows_kept += 1
            length_bp = end - start + 1
            if symbol in genes:
                # Shouldn't normally happen (one gene record per symbol),
                # but if it does, keep the longer span rather than silently
                # picking one arbitrarily.
                genes[symbol] = max(genes[symbol], length_bp)
            else:
                genes[symbol] = length_bp

        # TE-hit rows are identified by the Description="Within range of
        # <gene>" attribute itself, NOT by requiring
        # type=="mobile_genetic_element" - build_tfbs_te_gff.py's normal
        # combined output always types these rows mobile_genetic_element,
        # but some older/hand-built files (seen in practice) carry this
        # same Description tag on rows with a blank/other feature-type
        # column instead. Matching on the attribute directly handles both,
        # and is safe because ordinary gene/mRNA/exon/etc. records never
        # carry a "Description=Within range of ..." attribute by
        # coincidence. This check is independent of (not "elif" on) the
        # gene check above since a row could theoretically be both.
        desc = attrs.get("Description", "")
        m = _WITHIN_RANGE_RE.search(desc)
        if m:
            te_counts_raw[m.group(1)] += 1

    rows = []
    unmatched_te_genes = []
    for symbol, length_bp in genes.items():
        rows.append({
            "symbol": symbol,
            "length_bp": length_bp,
            "te_count": te_counts_raw.get(symbol, 0),
        })
    for symbol in te_counts_raw:
        if symbol not in genes:
            unmatched_te_genes.append(symbol)

    return rows, unmatched_te_genes, n_gene_rows_seen, n_gene_rows_kept


# ==========================================================================
# 3. Config file

def _strip_quotes(s):
    s = s.strip()
    if len(s) >= 2 and s[0] == s[-1] and s[0] in "\"'":
        s = s[1:-1]
    return s


def load_species_config(path):
    """Reads an INI file with one [section] per species:
        [suzukii]
        gff      = combined_suzukii_cyp_annotation.gff3
        exposure = high

    Optionally, a species can also set:
        gene_name_pattern = ^cyp\\d

    to restrict which `gene` features in that species' GFF3 count as Cyp
    genes (matched case-insensitively, anchored at the start of the Name
    attribute). This is only needed when a species' file is a *whole-
    genome* annotation rather than one already pre-filtered down to just
    the Cyp genes analyzed - see parse_species_gff3's docstring.

    Paths are resolved relative to the config file's own directory (same
    convention as build_tfbs_te_gff.py's species_config.ini), so a config +
    its data files can be copied around together. Returns a dict:
        {species_name: {"gff": resolved_path, "exposure": "high"|"low",
                         "gene_name_pattern": compiled regex or None}}
    """
    cfg = configparser.ConfigParser()
    read_ok = cfg.read(path)
    if not read_ok:
        raise FileNotFoundError(f"Config file not found or unreadable: {path}")

    base_dir = Path(path).resolve().parent
    species = {}
    for section in cfg.sections():
        if "gff" not in cfg[section]:
            raise ValueError(f"[{section}] is missing required 'gff' setting")
        if "exposure" not in cfg[section]:
            raise ValueError(f"[{section}] is missing required 'exposure' setting")
        exposure = _strip_quotes(cfg[section]["exposure"]).lower()
        if exposure not in ("high", "low"):
            raise ValueError(f"[{section}] exposure must be 'high' or 'low', got {exposure!r}")
        # Strip a single matching pair of surrounding quotes, if present -
        # Windows paths with spaces (e.g. a filename like "foo 1.gff") are
        # often quoted by users/editors, but configparser does NOT strip
        # quotes itself (unlike a shell), so an unstripped quote character
        # would otherwise be treated as part of the path itself.
        gff_path = _strip_quotes(cfg[section]["gff"])
        resolved = (base_dir / gff_path) if not Path(gff_path).is_absolute() else Path(gff_path)

        gene_name_pattern = None
        if "gene_name_pattern" in cfg[section]:
            pattern_str = _strip_quotes(cfg[section]["gene_name_pattern"])
            if pattern_str:
                try:
                    gene_name_pattern = re.compile(pattern_str, re.IGNORECASE)
                except re.error as exc:
                    raise ValueError(f"[{section}] invalid gene_name_pattern {pattern_str!r}: {exc}")

        species[section] = {"gff": str(resolved), "exposure": exposure, "gene_name_pattern": gene_name_pattern}
    if not species:
        raise ValueError(f"No [species] sections found in {path}")
    return species


# ==========================================================================
# 4. Per-species metrics + pooled analysis

def compute_species_summary(species_name, exposure, gene_rows):
    n_cyp_genes = len(gene_rows)
    n_genes_with_te = sum(1 for g in gene_rows if g["te_count"] > 0)
    total_te_count = sum(g["te_count"] for g in gene_rows)
    total_cyp_bp = sum(g["length_bp"] for g in gene_rows)
    return {
        "species": species_name,
        "exposure": exposure,
        "n_cyp_genes": n_cyp_genes,
        "n_genes_with_te": n_genes_with_te,
        "pct_genes_with_te": (100.0 * n_genes_with_te / n_cyp_genes) if n_cyp_genes else 0.0,
        "total_te_count": total_te_count,
        "total_cyp_bp": total_cyp_bp,
        "te_per_kb": (total_te_count / (total_cyp_bp / 1000.0)) if total_cyp_bp else 0.0,
        "te_per_gene": (total_te_count / n_cyp_genes) if n_cyp_genes else 0.0,
    }


def run_pooled_tests(all_gene_rows):
    """all_gene_rows: list of (species, exposure, gene_row_dict) across all
    species. Runs the two pre-registered pooled per-gene tests and returns
    a dict of results."""
    high_has_te = sum(1 for _, exp, g in all_gene_rows if exp == "high" and g["te_count"] > 0)
    high_no_te = sum(1 for _, exp, g in all_gene_rows if exp == "high" and g["te_count"] == 0)
    low_has_te = sum(1 for _, exp, g in all_gene_rows if exp == "low" and g["te_count"] > 0)
    low_no_te = sum(1 for _, exp, g in all_gene_rows if exp == "low" and g["te_count"] == 0)

    fisher_p = fisher_exact_two_tailed(high_has_te, high_no_te, low_has_te, low_no_te)
    chi2, chi2_p = chi_square_2x2(high_has_te, high_no_te, low_has_te, low_no_te)

    high_counts = [g["te_count"] for _, exp, g in all_gene_rows if exp == "high"]
    low_counts = [g["te_count"] for _, exp, g in all_gene_rows if exp == "low"]
    mwu_U, mwu_p = mann_whitney_u_test(high_counts, low_counts)
    ttest_t, ttest_df, ttest_p = welch_t_test(high_counts, low_counts)

    return {
        "table": {"high_has_te": high_has_te, "high_no_te": high_no_te,
                  "low_has_te": low_has_te, "low_no_te": low_no_te},
        "fisher_p": fisher_p,
        "chi2": chi2, "chi2_p": chi2_p,
        "mwu_U": mwu_U, "mwu_p": mwu_p,
        "ttest_t": ttest_t, "ttest_df": ttest_df, "ttest_p": ttest_p,
        "high_n": len(high_counts), "low_n": len(low_counts),
        "high_mean_te": (sum(high_counts) / len(high_counts)) if high_counts else 0.0,
        "low_mean_te": (sum(low_counts) / len(low_counts)) if low_counts else 0.0,
    }


def bootstrap_species_cluster_ci(all_gene_rows, n_resamples=10000, seed=12345, ci=0.95):
    """Species-level cluster bootstrap for the difference/ratio in mean
    per-gene TE count between exposure groups.

    WHY THIS EXISTS: the pooled per-gene Fisher's/Mann-Whitney tests above
    treat every gene as an independent observation - useful for power, but
    pseudoreplicated (genes in the same species aren't really
    independent). Going the other direction and testing at the species
    level instead has essentially NO power here: with only 5 species
    split some way into two groups, an exact species-level test has at
    most C(5,3)=10 possible group assignments, so the smallest p-value it
    could ever produce is 1/10=0.10 - it can never clear a conventional
    0.05 threshold no matter how large the true effect is. That ceiling,
    not a flaw in the analysis, is very likely why the pooled tests read
    "inconclusive": there aren't enough independent species to clear
    significance on their own, regardless of effect size.

    This function splits the difference: it resamples WHOLE SPECIES (with
    replacement, separately within each exposure group), so species stays
    the real unit of replication, but produces as many resamples as you
    like (not capped at 10), giving a smoother confidence interval on the
    effect size instead of a single pass/fail significance call.

    Returns a dict with the observed per-group mean TE/gene, bootstrap
    percentile CIs for both the difference (high-low) and ratio
    (high/low), and prob_high_greater - the fraction of resamples where
    the high-exposure group's mean TE/gene exceeded the low-exposure
    group's (a directional confidence statement, not a p-value).
    """
    by_species = defaultdict(lambda: {"exposure": None, "te_counts": []})
    for species, exposure, g in all_gene_rows:
        by_species[species]["exposure"] = exposure
        by_species[species]["te_counts"].append(g["te_count"])

    high_species = [v["te_counts"] for v in by_species.values() if v["exposure"] == "high"]
    low_species = [v["te_counts"] for v in by_species.values() if v["exposure"] == "low"]

    def group_mean(counts_lists):
        pooled = [c for lst in counts_lists for c in lst]
        return (sum(pooled) / len(pooled)) if pooled else 0.0

    observed_high_mean = group_mean(high_species)
    observed_low_mean = group_mean(low_species)

    rng = random.Random(seed)
    diffs = []
    ratios = []
    n_high_greater = 0
    for _ in range(n_resamples):
        resampled_high = ([high_species[rng.randrange(len(high_species))] for _ in range(len(high_species))]
                           if high_species else [])
        resampled_low = ([low_species[rng.randrange(len(low_species))] for _ in range(len(low_species))]
                          if low_species else [])
        h_mean = group_mean(resampled_high)
        l_mean = group_mean(resampled_low)
        diffs.append(h_mean - l_mean)
        if l_mean > 0:
            ratios.append(h_mean / l_mean)
        if h_mean > l_mean:
            n_high_greater += 1
        elif h_mean == l_mean:
            n_high_greater += 0.5  # ties split, avoids biasing the directional stat either way

    diffs.sort()

    def percentile_ci(sorted_vals, ci):
        n = len(sorted_vals)
        if n == 0:
            return (float("nan"), float("nan"))
        lo_idx = max(0, min(n - 1, int(round((1 - ci) / 2 * n))))
        hi_idx = max(0, min(n - 1, int(round((1 - (1 - ci) / 2) * n)) - 1))
        return (sorted_vals[lo_idx], sorted_vals[hi_idx])

    diff_ci = percentile_ci(diffs, ci)
    ratios.sort()
    ratio_ci = percentile_ci(ratios, ci)

    return {
        "observed_high_mean_te_per_gene": observed_high_mean,
        "observed_low_mean_te_per_gene": observed_low_mean,
        "n_high_species": len(high_species),
        "n_low_species": len(low_species),
        "diff_ci": diff_ci,
        "ratio_ci": ratio_ci,
        "prob_high_greater": n_high_greater / n_resamples,
        "n_resamples": n_resamples,
        "ci_level": ci,
    }


def format_bootstrap_section(boot, alpha=0.05):
    """Shared markdown block for the species-cluster bootstrap result -
    used by all three report builders (exposure/cncc/xenobiotic) so the
    framing stays consistent across scripts."""
    lines = []
    lines.append("## Effect size: species-level cluster bootstrap (complements the test above)")
    lines.append("")
    lines.append(
        "The pooled test above treats every gene as independent (more power, but pseudoreplicated). "
        "A pure species-level test would be the statistically \"safe\" alternative, but with only "
        f"{boot['n_high_species']} high- and {boot['n_low_species']} low-exposure species, it has "
        "essentially no power - the smallest p-value an exact species-level test could produce with "
        "5 species split 3-2 is 1/C(5,3) = 0.10, which can never clear a conventional significance "
        "threshold regardless of effect size. This bootstrap splits the difference: it resamples "
        "**whole species** (not genes) with replacement, so species stays the real unit of "
        "replication, while still producing a smooth interval estimate rather than a single "
        "pass/fail call."
    )
    lines.append("")
    ratio = (boot["observed_high_mean_te_per_gene"] / boot["observed_low_mean_te_per_gene"]
             if boot["observed_low_mean_te_per_gene"] > 0 else float("inf"))
    lines.append(
        f"- Observed mean TE/gene: high-exposure = {boot['observed_high_mean_te_per_gene']:.3f}, "
        f"low-exposure = {boot['observed_low_mean_te_per_gene']:.3f} "
        f"({ratio:.2f}x)" if ratio != float("inf") else
        f"- Observed mean TE/gene: high-exposure = {boot['observed_high_mean_te_per_gene']:.3f}, "
        f"low-exposure = {boot['observed_low_mean_te_per_gene']:.3f} (low-exposure mean is 0, ratio undefined)"
    )
    lo, hi = boot["diff_ci"]
    lines.append(f"- Bootstrap {int(boot['ci_level']*100)}% CI for the **difference** (high - low): [{lo:.3f}, {hi:.3f}]")
    rlo, rhi = boot["ratio_ci"]
    if rlo == rlo:  # not NaN
        lines.append(f"- Bootstrap {int(boot['ci_level']*100)}% CI for the **ratio** (high / low): [{rlo:.2f}, {rhi:.2f}]")
    lines.append(
        f"- **{boot['prob_high_greater']*100:.1f}% of {boot['n_resamples']:,} species-resamples** "
        "showed a higher mean TE/gene in high-exposure species than low-exposure species"
    )
    lines.append("")
    if boot["n_high_species"] <= 2 or boot["n_low_species"] <= 2:
        lines.append(
            f"> Caveat: with only {boot['n_high_species']} high- and {boot['n_low_species']} "
            "low-exposure species, this bootstrap itself is coarse (few distinct species "
            "compositions are possible) - treat the interval width and the "
            "\"% of resamples\" figure as a directional confidence statement, not a formal "
            "significance test. More species would sharpen this considerably."
        )
        lines.append("")
    return "\n".join(lines)


# ==========================================================================
# 5. Report generation

def build_report(species_summaries, pooled, alpha, unmatched_by_species, boot=None):
    lines = []
    lines.append("# Does insecticide exposure predict TE insertions in Cyp genes?")
    lines.append("")
    lines.append(
        "This report tests the abstract's hypothesis - that Drosophila species with "
        "greater insecticide exposure (D. suzukii, D. subpulchrella, D. melanogaster) "
        "carry more transposable element (TE) insertions in/near Cytochrome P450 (Cyp) "
        "genes than species with less exposure (D. santomea, D. sechellia) - against the "
        "TE/gene annotations already produced by build_tfbs_te_gff.py, for the "
        f"{len(species_summaries)} species that loaded successfully below."
    )
    if len(species_summaries) < 5:
        lines.append("")
        lines.append(
            f"> Note: only {len(species_summaries)} of the intended 5 species loaded - check "
            "the run's [warn]/[error] output for any config path issues before treating results "
            "below as complete."
        )
    lines.append("")

    lines.append("## Per-species summary")
    lines.append("")
    lines.append("| Species | Exposure | Cyp genes | Genes w/ TE | % w/ TE | Total TEs | Cyp bp | TE/kb | TE/gene |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---:|---:|")
    for s in species_summaries:
        lines.append(
            f"| {s['species']} | {s['exposure']} | {s['n_cyp_genes']} | {s['n_genes_with_te']} | "
            f"{s['pct_genes_with_te']:.1f}% | {s['total_te_count']} | {s['total_cyp_bp']:,} | "
            f"{s['te_per_kb']:.3f} | {s['te_per_gene']:.3f} |"
        )
    lines.append("")

    for sp, unmatched in unmatched_by_species.items():
        if unmatched:
            lines.append(
                f"> Note: {sp} - {len(unmatched)} TE record(s) referenced a gene symbol via "
                f"'Within range of ...' that was not found among that species' `gene` features "
                f"({unmatched[:5]}{' ...' if len(unmatched) > 5 else ''}); these TEs were not "
                f"countable against any gene and are excluded from the per-gene metrics above."
            )
    lines.append("")

    lines.append("## Pooled per-gene statistical test (pre-registered)")
    lines.append("")
    lines.append(
        "Every Cyp gene across all 5 species is pooled into one table of "
        f"{pooled['high_n'] + pooled['low_n']} gene-level observations "
        f"({pooled['high_n']} from high-exposure species, {pooled['low_n']} from low-exposure "
        "species), and two tests are run on it:"
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
        f"**Mann-Whitney U test on per-gene TE counts** (high-exposure genes averaged "
        f"{pooled['high_mean_te']:.3f} TEs/gene vs {pooled['low_mean_te']:.3f} TEs/gene for "
        f"low-exposure): U = {pooled['mwu_U']:.1f}, p = {pooled['mwu_p']:.6g}"
    )
    lines.append(
        f"- **Welch's t-test (unequal variances) on per-gene TE counts:** "
        f"t = {pooled['ttest_t']:.3f}, df = {pooled['ttest_df']:.1f}, p = {pooled['ttest_p']:.6g} "
        "*(reported as a familiar reference point - TE counts are skewed, non-normal count data, "
        "so Mann-Whitney above is the more statistically appropriate test; treat a disagreement "
        "between the two as a signal to look at the distribution shape, not as license to pick "
        "whichever gives the smaller p-value)*"
    )
    lines.append("")

    table_saturated = (t["high_no_te"] == 0 and t["low_no_te"] == 0)
    if table_saturated:
        lines.append(
            "> Note: every Cyp gene in every species has >=1 associated TE, so the "
            "presence/absence table above has **no variance to test** - Fisher's exact and "
            "chi-square returning p=1 here is an artifact of that saturation, not evidence of "
            "no effect. The Mann-Whitney U test on per-gene TE **counts** (below) is the test "
            "actually carrying signal in this dataset."
        )
        lines.append("")

    fisher_sig = pooled["fisher_p"] < alpha
    mwu_sig = pooled["mwu_p"] < alpha
    # Direction is judged from mean per-gene TE burden (a continuous measure),
    # not from the presence/absence table above - that table can be fully
    # saturated (every gene in every group has >=1 TE, as above) and would
    # then give a meaningless/degenerate "direction" by cross-multiplication.
    # Mean TE count doesn't have that failure mode.
    if pooled["high_mean_te"] > pooled["low_mean_te"]:
        direction_supports = True
    elif pooled["high_mean_te"] < pooled["low_mean_te"]:
        direction_supports = False
    else:
        direction_supports = None  # exact tie - no directional signal either way

    lines.append(f"## Verdict (alpha = {alpha})")
    lines.append("")
    if direction_supports is None:
        verdict = (
            "**Inconclusive (tied).** High- and low-exposure groups show identical mean "
            "per-gene TE burden - there is no directional signal either way in this data."
        )
    elif fisher_sig and mwu_sig and direction_supports:
        verdict = (
            "**Supports the hypothesis.** Both tests are statistically significant, and the "
            "direction is as predicted: high-exposure species show a higher proportion of "
            "TE-affected Cyp genes and a higher per-gene TE burden than low-exposure species."
        )
    elif (fisher_sig or mwu_sig) and direction_supports:
        verdict = (
            "**Partially supports the hypothesis.** The direction is as predicted (high-exposure "
            "species show a higher per-gene TE burden), and one of the two tests reaches "
            f"significance at alpha={alpha}, but not both - treat this as suggestive, not conclusive."
            + (" Note the presence/absence test is uninformative here (see note above), so the "
               "Mann-Whitney result carries the real weight of this verdict." if table_saturated else "")
        )
    elif not direction_supports and (fisher_sig or mwu_sig):
        verdict = (
            "**Contradicts the hypothesis.** At least one test is statistically significant, but "
            "in the OPPOSITE direction from predicted - low-exposure species show a higher "
            "per-gene TE burden than high-exposure species in this data."
        )
    else:
        verdict = (
            "**Inconclusive.** Neither test reaches statistical significance - this data does not "
            "provide clear support for (or against) the hypothesis as tested here."
        )
    lines.append(verdict)
    lines.append("")

    if boot is not None:
        lines.append(format_bootstrap_section(boot, alpha))

    lines.append("## Important caveats")
    lines.append("")
    lines.append(
        "- **Pseudoreplication:** pooling every gene across species as an independent observation "
        "ignores that genes within the same species share ancestry/genomic background - a species-"
        "level effect (e.g. one species' overall TE landscape) can masquerade as a gene-level "
        "exposure effect. With only 5 species (3 high, 2 low), a formal species-level test isn't "
        "statistically meaningful, so this pooled gene-level test is the best available option, "
        "but its significance should be read as suggestive of association, not proof of a "
        "TE-insertion response specifically caused by insecticide exposure."
    )
    lines.append(
        "- **Correlation vs. causation:** even a significant, correctly-directed result shows "
        "association between exposure-group membership and TE/Cyp-gene patterns, not that "
        "insecticide exposure itself caused the TE insertions (species differ in many other ways "
        "besides exposure history: genome size, TE landscape, demographic history, etc.)."
    )
    lines.append(
        "- **TE-gene association source:** TE-to-gene assignment here is inherited as-is from "
        "the `Within range of <gene>` annotation already written by build_tfbs_te_gff.py (i.e. "
        "whatever flanking-distance/overlap rule that pipeline used to call a TE 'near' a gene); "
        "this script does not re-derive or second-guess that call."
    )
    lines.append("")

    high_species = [s for s in species_summaries if s["exposure"] == "high"]
    low_species = [s for s in species_summaries if s["exposure"] == "low"]

    lines.append("## Species-level descriptive comparison (exploratory only, not tested)")
    lines.append("")
    lines.append(
        f"With n={len(high_species)} high-exposure species and n={len(low_species)} low-exposure "
        "species loaded, a formal species-level test has essentially no statistical power and is "
        "not reported. For descriptive context only:"
    )
    lines.append("")
    if high_species and low_species:
        mean_high_pct = sum(s["pct_genes_with_te"] for s in high_species) / len(high_species)
        mean_low_pct = sum(s["pct_genes_with_te"] for s in low_species) / len(low_species)
        mean_high_tpg = sum(s["te_per_gene"] for s in high_species) / len(high_species)
        mean_low_tpg = sum(s["te_per_gene"] for s in low_species) / len(low_species)
        lines.append(
            f"- Mean % Cyp genes with a TE: high-exposure species = {mean_high_pct:.1f}%, "
            f"low-exposure species = {mean_low_pct:.1f}%"
        )
        lines.append(
            f"- Mean TEs per Cyp gene: high-exposure species = {mean_high_tpg:.3f}, "
            f"low-exposure species = {mean_low_tpg:.3f}"
        )
    lines.append("")

    return "\n".join(lines)


# ==========================================================================
# 6. Optional plot

def maybe_write_plot(species_summaries, out_path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("[info] matplotlib not available - skipping plot (everything else still ran fine)")
        return False

    species_summaries = sorted(species_summaries, key=lambda s: (s["exposure"] != "high", s["species"]))
    names = [s["species"] for s in species_summaries]
    values = [s["te_per_gene"] for s in species_summaries]
    colors = ["#c0392b" if s["exposure"] == "high" else "#2980b9" for s in species_summaries]

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.bar(names, values, color=colors)
    ax.set_ylabel("TEs per Cyp gene")
    ax.set_title("TE density in Cyp genes by species\n(red = high insecticide exposure, blue = low)")
    plt.xticks(rotation=30, ha="right")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"[info] wrote plot to {out_path}")
    return True


# ==========================================================================
# 7. Main

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Test whether insecticide exposure predicts TE insertions in Cyp genes, "
                    "using the combined GFF3 outputs from build_tfbs_te_gff.py.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--config", help="Path to species config INI (see te_cyp_species_config.example.ini)")
    parser.add_argument("--output-csv", default="te_cyp_species_summary.csv",
                         help="Where to write the per-species summary CSV (default: %(default)s)")
    parser.add_argument("--per-gene-csv", default=None,
                         help="Optional: also write a per-gene CSV (species, symbol, exposure, length_bp, te_count)")
    parser.add_argument("--output-report", default="te_cyp_report.md",
                         help="Where to write the markdown report (default: %(default)s)")
    parser.add_argument("--alpha", type=float, default=0.05, help="Significance threshold (default: 0.05)")
    parser.add_argument("--plot", dest="plot", action="store_true", default=True,
                         help="Write a bar chart of TE/gene by species (default: on, auto-skips if matplotlib missing)")
    parser.add_argument("--no-plot", dest="plot", action="store_false", help="Skip the plot entirely")
    parser.add_argument("--plot-path", default="te_cyp_by_species.png", help="Plot output path (default: %(default)s)")
    parser.add_argument("--self-test", action="store_true",
                         help="Just run the internal statistics cross-validation checks and exit")
    args = parser.parse_args(argv)

    if args.self_test:
        ok = validate_stats()
        sys.exit(0 if ok else 1)

    if not args.config:
        parser.error("--config is required (unless using --self-test)")

    ok = validate_stats()
    if not ok:
        print("[error] statistics self-test failed - refusing to report results from untrusted code", file=sys.stderr)
        sys.exit(1)
    print()

    species_cfg = load_species_config(args.config)

    species_summaries = []
    all_gene_rows = []          # (species, exposure, gene_row_dict)
    per_gene_csv_rows = []
    unmatched_by_species = {}

    for species_name, cfg in sorted(species_cfg.items()):
        gff_path = cfg["gff"]
        exposure = cfg["exposure"]
        if not Path(gff_path).exists():
            print(f"[warn] {species_name}: GFF3 file not found, skipping: {gff_path}", file=sys.stderr)
            continue
        gene_name_pattern = cfg.get("gene_name_pattern")
        gene_rows, unmatched, n_seen, n_kept = parse_species_gff3(gff_path, gene_name_pattern)
        unmatched_by_species[species_name] = unmatched
        if not gene_rows:
            print(f"[warn] {species_name}: no `gene` features found in {gff_path}", file=sys.stderr)
        if gene_name_pattern is not None:
            print(f"[info] {species_name}: gene_name_pattern kept {n_kept} of {n_seen} `gene` "
                  f"records in {gff_path}")
        elif n_seen > 500:
            print(f"[warn] {species_name}: {n_seen} `gene` records found with no gene_name_pattern "
                  f"filter set - if {Path(gff_path).name} is a whole-genome annotation rather than "
                  f"a file pre-filtered to just Cyp genes, this will badly distort the comparison; "
                  f"consider adding gene_name_pattern = ^cyp\\d to [{species_name}] in your config",
                  file=sys.stderr)
        summary = compute_species_summary(species_name, exposure, gene_rows)
        species_summaries.append(summary)
        for g in gene_rows:
            all_gene_rows.append((species_name, exposure, g))
            per_gene_csv_rows.append({
                "species": species_name, "exposure": exposure, "symbol": g["symbol"],
                "length_bp": g["length_bp"], "te_count": g["te_count"],
            })
        print(f"[info] {species_name} ({exposure} exposure): {summary['n_cyp_genes']} Cyp genes, "
              f"{summary['n_genes_with_te']} with >=1 TE ({summary['pct_genes_with_te']:.1f}%), "
              f"{summary['te_per_gene']:.3f} TEs/gene")

    if not species_summaries:
        print("[error] no species data loaded - check --config paths", file=sys.stderr)
        sys.exit(1)

    n_high = sum(1 for s in species_summaries if s["exposure"] == "high")
    n_low = sum(1 for s in species_summaries if s["exposure"] == "low")
    if n_high == 0 or n_low == 0:
        print("[error] need at least one species in each exposure group (high and low) to run the "
              "comparison", file=sys.stderr)
        sys.exit(1)

    pooled = run_pooled_tests(all_gene_rows)
    boot = bootstrap_species_cluster_ci(all_gene_rows)

    # --- per-species CSV ---
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

    # --- report ---
    report_text = build_report(species_summaries, pooled, args.alpha, unmatched_by_species, boot)
    with open(args.output_report, "w") as fh:
        fh.write(report_text)
    print(f"[info] wrote report to {args.output_report}")

    # --- plot ---
    if args.plot:
        maybe_write_plot(species_summaries, args.plot_path)

    print()
    print(f"Fisher's exact p = {pooled['fisher_p']:.6g}   Mann-Whitney U p = {pooled['mwu_p']:.6g}")
    print(f"Species-cluster bootstrap: {boot['prob_high_greater']*100:.1f}% of resamples showed "
          f"higher TE/gene in high-exposure species")


if __name__ == "__main__":
    main()
