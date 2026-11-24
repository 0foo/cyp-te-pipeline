# NEW (not a lab file, 2026-11-23): builds summary/pipeline-results.xlsx (uploaded as the Google Sheet
# "Cyp-TE pipeline results"). It reads only this pipeline's outputs: summary/per-species.tsv (stages 1-3,
# from make-summary.py) and the stage 4 outputs in 4-combine-and-compare/output/compare/ (repo scripts
# run unchanged by run-compare.sh).
# TABS: About; Per species (the stage 1-3 summary joined with each compare script's per-species
#   summary); Statistics (one row per comparison, giving the Fisher, chi-square, Mann-Whitney and Welch
#   results, the bootstrap CIs and the verdict); Cyp genes (stage 4) and Xenobiotic genes (stage 4) (the
#   scripts' per-gene CSVs). Each data tab has a "<tab> - columns" tab defining every column.
# The statistics are copied from the scripts' own markdown reports by regular expressions, not
#   recomputed. If a report's wording changes, a field comes out blank instead of raising an error.
# Usage: nix-shell -p 'python3.withPackages(p: [p.openpyxl])' --run 'python3 -B make-sheet.py'
import csv, re
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment
from openpyxl.utils import get_column_letter

A = "/bioinformatics/atallah-repeats"
C = f"{A}/4-combine-and-compare/output/compare"
TESTS = ["exposure", "cncc", "xenobiotic"]
NAMES = {"exposure": "All Cyp genes", "cncc": "CncC-proximal Cyp genes", "xenobiotic": "Xenobiotic-list Cyp genes"}

def num(v):
    try:
        f = float(v); return int(f) if f.is_integer() and "." not in v else round(f, 4)
    except (ValueError, TypeError):
        return v

def tsv(p, d="\t"): return list(csv.DictReader(open(p), delimiter=d))

def sheet(wb, title, rows, cols=None, widths=None):
    ws = wb.create_sheet(title)
    cols = cols or list(rows[0].keys())
    ws.append(cols)
    for r in rows: ws.append([num(r.get(c, "")) for c in cols])
    for c in ws[1]: c.font = Font(bold=True)
    ws.freeze_panes = "B2"; ws.auto_filter.ref = ws.dimensions
    for i, c in enumerate(cols, 1):
        ws.column_dimensions[get_column_letter(i)].width = (widths or {}).get(c, max(10, min(40, len(c) + 2)))
    return ws

# ---- statistics, parsed from the reports
def stats(t):
    r = open(f"{C}/{t}.report.md").read()
    g = lambda pat: (m.group(1) if (m := re.search(pat, r, re.S)) else "")
    tab = re.findall(r"\| (High|Low) exposure \| (\d+) \| (\d+) \|", r)
    return {
        "comparison": NAMES[t], "script": f"compare_te_cyp_{t}.py",
        "high_species": g(r"only (\d+) high- and \d+ low-exposure species"),
        "low_species": g(r"only \d+ high- and (\d+) low-exposure species"),
        "high_genes": g(r"\((\d+) (?:from )?high-exposure"), "low_genes": g(r"(\d+) (?:from )?low-exposure(?: species)?\)"),
        "high_genes_with_te / without": f"{tab[0][1]} / {tab[0][2]}" if tab else "",
        "low_genes_with_te / without": f"{tab[1][1]} / {tab[1][2]}" if len(tab) > 1 else "",
        "fisher_p": g(r"Fisher's exact test \(two-tailed\):\*\* p = (\S+)"),
        "chi2_yates": g(r"chi2 = (\S+),"), "chi2_p": g(r"chi2 = \S+, p = (\S+)"),
        "high_mean_te_per_gene": g(r"averaged ([\d.]+) TEs/gene vs"), "low_mean_te_per_gene": g(r"vs ([\d.]+) TEs/gene for low"),
        "mannwhitney_U": g(r"U = ([\d.]+), p"), "mannwhitney_p": g(r"U = [\d.]+, p = (\S+)"),
        "welch_t": g(r"t = (\S+), df"), "welch_df": g(r"df = (\S+), p"), "welch_p": g(r"df = \S+, p = (\S+)"),
        "bootstrap_ratio_high_low": g(r"low-exposure = [\d.]+ \(([\d.]+)x\)"),
        "bootstrap_diff_95ci": g(r"\*\*difference\*\* \(high - low\): \[([^\]]+)\]"),
        "bootstrap_ratio_95ci": g(r"\*\*ratio\*\* \(high / low\): \[([^\]]+)\]"),
        "pct_resamples_high_gt_low": g(r"\*\*([\d.]+%) of [\d,]+ species-resamples"),
        "verdict": g(r"## Verdict[^\n]*\n\n\*\*([^*]+)\*\*"),
    }
st = [stats(t) for t in TESTS]

# ---- per species: stages 1-3 summary + stage 4 per-species summaries
sp = tsv(f"{A}/summary/per-species.tsv")
s4 = {t: {r["species"]: r for r in tsv(f"{C}/{t}.species-summary.csv", ",")} for t in TESTS}
S4COLS = ["n_cyp_genes", "n_genes_with_te", "total_te_count", "total_cyp_bp", "te_per_kb", "te_per_gene"]
for r in sp:
    for t in TESTS:
        x = s4[t].get(r["species"], {})
        if t == "cncc": r["s4_cncc_n_tfbs_total"] = x.get("n_tfbs_total", "")
        for c in S4COLS: r[f"s4_{t}_{c}"] = x.get(c, "")


# ---- column dictionaries: one "<tab> - columns" tab per data tab
W = "stage 3 window = Cyp gene +/-3 kb (Locate_TE.py: the repeat must lie fully inside)"
S4 = ("Counted by the compare script from the species' combined GFF3: a 'TE' is any mobile_genetic_element row whose "
      "Description says 'Within range of <gene>' (every stage 3 RepeatMasker hit, simple repeats included).")
S4D = {
    "n_cyp_genes": "Number of Cyp gene records (unique Name=) the script counted in this species.",
    "n_genes_with_te": "How many of those genes have >=1 TE.",
    "total_te_count": "Total TEs across those genes.",
    "total_cyp_bp": "Summed length (bp) of those genes (gene start to end).",
    "te_per_kb": "total_te_count / (total_cyp_bp / 1000).",
    "te_per_gene": "total_te_count / n_cyp_genes. This is the per-species value the tests compare.",
}
SUBSET = {"exposure": "all Cyp genes", "cncc": "Cyp genes near a predicted CncC:Maf-S site only",
          "xenobiotic": "the 10 xenobiotic-list genes only"}
PER_SPECIES = {
    "species": ("ID", "Species label used in this pipeline (annotation/GFF label, upper case)."),
    "rm_sample": ("ID", "Genome sample name used for RepeatModeler/RepeatMasker and the genome FASTA (.nanopore etc. = which assembly)."),
    "sheet_label": ("Ecology sheet", "Row name in the ecology Google Sheet this species was matched to."),
    "exposure_group": ("Ecology sheet -> make-exposure.py", "high = 'Documented exposure...' or 'Likely'; low = 'Unlikely'/'Very unlikely'; none = Possible/Possible to likely/Unknown (left out of all tests)."),
    "insecticide_exposure": ("Ecology sheet", "The sheet's own 'Insecticide exposure' call (sheet 3, column E) that exposure_group is based on."),
    "exposure_setting": ("Ecology sheet", "Where exposure would happen (orchards, households, wild habitat...)."),
    "resistance_evidence": ("Ecology sheet", "Published field/lab insecticide resistance, if any."),
    "confounder": ("make-exposure.py", "'host-toxin confounder' for low species whose detox genes track plant/fungal toxins (cactophiles, wild-fungus breeders, D. sechellia); blank otherwise."),
    "genome_bp": ("Stage 1 (.rm.tbl)", "Total assembly length, bp."),
    "genome_seqs": ("Stage 1 (.rm.tbl)", "Number of sequences (contigs/scaffolds) in the assembly."),
    "pct_masked": ("Stage 1 (.rm.tbl)", "% of the genome masked by RepeatMasker (all repeat types)."),
    "pct_LINE": ("Stage 1 (.rm.tbl)", "% of genome in LINEs."),
    "pct_LTR": ("Stage 1 (.rm.tbl)", "% of genome in LTR retrotransposons."),
    "pct_DNA": ("Stage 1 (.rm.tbl)", "% of genome in DNA transposons."),
    "pct_RC": ("Stage 1 (.rm.tbl)", "% of genome in rolling-circle (Helitron) elements."),
    "pct_unclassified": ("Stage 1 (.rm.tbl)", "% of genome in repeats RepeatModeler could not classify."),
    "pct_simple": ("Stage 1 (.rm.tbl)", "% of genome in simple repeats (microsatellites)."),
    "repeat_families": ("Stage 1 (families.fa)", "Number of repeat families in the species' RepeatModeler library."),
    "gene_frac_on_rm_seqs": ("Stage 3 check (make-pairs.py)", "Fraction of annotated genes whose sequence is in the RepeatMasker output: checks the annotation and the repeat run use the same assembly (>=0.90 to run)."),
    "cyp_genes_in_annotation": ("Stage 3 check (make-pairs.py)", "Gene records in the renamed annotation whose attributes contain a Cyp name from Reg_Gene_Full.txt (substring match, as the lab scripts do)."),
    "cyp_loci": ("Stage 3 (repeatOpp filtered.gff)", "Distinct Cyp loci (unique gene ID) that repeatOpp.py kept."),
    "cyp_loci_rows_raw": ("Stage 3 (repeatOpp filtered.gff)", "Rows in filtered.gff before de-duplication (repeatOpp emits duplicate rows, a preserved lab defect)."),
    "te_table_rows_raw": ("Stage 3 (GenesAffectedByTEs.txt)", "Rows in the stage 3 gene x repeat table, duplicates included."),
    "unique_hits": ("Stage 3", "Distinct (Cyp locus, RepeatMasker row) pairs: every repeat hit in a " + W + "."),
    "te_hits": ("Stage 3", "unique_hits that are transposable elements (excludes Simple_repeat, Low_complexity, Satellite, rRNA/tRNA/snRNA/scRNA/srpRNA, ARTEFACT)."),
    "simple_lowcomplex_hits": ("Stage 3", "unique_hits that are not TEs (the excluded classes above)."),
    "pct_hits_te": ("Stage 3", "te_hits / unique_hits x 100."),
    "cyp_loci_with_te": ("Stage 3", "Cyp loci with >=1 TE hit."),
    "pct_cyp_loci_with_te": ("Stage 3", "cyp_loci_with_te / cyp_loci x 100."),
    "cyp_loci_with_any_hit": ("Stage 3", "Cyp loci with >=1 hit of any kind (TE or simple)."),
    "te_per_cyp_locus": ("Stage 3", "te_hits / cyp_loci."),
    "te_per_10kb_window": ("Stage 3", "TE hits per 10 kb of window (window = locus length + 6 kb)."),
    "te_hits_LINE": ("Stage 3", "TE hits of class LINE."),
    "te_hits_LTR": ("Stage 3", "TE hits of class LTR."),
    "te_hits_DNA": ("Stage 3", "TE hits of class DNA (DNA transposons)."),
    "te_hits_RC": ("Stage 3", "TE hits of class RC (Helitrons)."),
    "te_hits_SINE": ("Stage 3", "TE hits of class SINE."),
    "te_hits_Unknown": ("Stage 3", "TE hits RepeatModeler left unclassified (Unknown)."),
    "te_hits_other": ("Stage 3", "TE hits in any other class (e.g. Retroposon)."),
    "xeno_targets_present": ("Stage 3", "How many of the 10 xenobiotic-list genes are among this species' Cyp loci (exact name match)."),
    "xeno_targets_with_te": ("Stage 3", "How many of those have >=1 TE hit."),
    "xeno_te_hits": ("Stage 3", "Total TE hits across the xenobiotic-list loci."),
    "s4_cncc_n_tfbs_total": ("Stage 4 (compare_te_cyp_cncc.py)", "All TF binding sites FIMO found in this species (JASPAR insect motifs + CncC motif, p<1e-4, in 1,000 bp up / 200 bp down of each TSS and around TEs). 0 would mean the motif scan failed."),
}
for t in TESTS:
    for c in S4COLS:
        PER_SPECIES[f"s4_{t}_{c}"] = (f"Stage 4 (compare_te_cyp_{t}.py), {SUBSET[t]}", S4D[c] + " Blank for 'none' species (not in the config).")
STATS = {
    "comparison": ("", "Which gene set was compared."),
    "script": ("", "Repo script that ran it (unchanged)."),
    "high_species": ("", "High-exposure species with data for this gene set."),
    "low_species": ("", "Low-exposure species with data for this gene set."),
    "high_genes": ("", "Gene-level observations from high species (genes pooled across species)."),
    "low_genes": ("", "Gene-level observations from low species."),
    "high_genes_with_te / without": ("", "2x2 table row: high-species genes with >=1 TE / with none."),
    "low_genes_with_te / without": ("", "2x2 table row: low-species genes with >=1 TE / with none."),
    "fisher_p": ("", "Fisher's exact test (two-tailed) on the 2x2 table. p = 1 here because no gene has zero TEs: no variance, not evidence of no effect."),
    "chi2_yates": ("", "Chi-square with Yates correction on the same table (cross-check)."),
    "chi2_p": ("", "p-value for chi2_yates."),
    "high_mean_te_per_gene": ("", "Mean TEs per gene, high-species genes pooled."),
    "low_mean_te_per_gene": ("", "Mean TEs per gene, low-species genes pooled."),
    "mannwhitney_U": ("", "Mann-Whitney U statistic on per-gene TE counts, high vs low (the scripts' main test)."),
    "mannwhitney_p": ("", "Two-tailed p for U. '0' means the p-value is too small for the script to represent (rounded to 0), not exactly zero."),
    "welch_t": ("", "Welch's t (unequal variances) on per-gene TE counts; negative = high < low. Reference only, counts are skewed."),
    "welch_df": ("", "Welch-Satterthwaite degrees of freedom."),
    "welch_p": ("", "Two-tailed p for welch_t."),
    "bootstrap_ratio_high_low": ("", "Observed mean TE/gene, high / low. Below 1 = high-exposure species have fewer."),
    "bootstrap_diff_95ci": ("", "95% CI of (high - low) mean TE/gene from 10,000 resamples of whole species (species, not genes, as the unit)."),
    "bootstrap_ratio_95ci": ("", "95% CI of the high/low ratio from the same bootstrap. Not containing 1 = a clear difference."),
    "pct_resamples_high_gt_low": ("", "Share of the 10,000 species resamples where high > low."),
    "verdict": ("", "The script's own verdict at alpha = 0.05 against the hypothesis (more TEs near Cyps in high-exposure species)."),
}
PER_GENE = {
    "species": ("", "Species label (high/low species only)."),
    "exposure": ("", "high or low (from exposure_group)."),
    "symbol": ("", "Cyp gene name (Name= in the combined GFF3, i.e. the D. melanogaster name assigned in stage 2)."),
    "length_bp": ("", "Gene length, bp (start to end; longest record if a name occurs twice)."),
    "te_count": ("", "TEs the script assigned to this gene. " + S4),
}

def doc_tab(wb, data_title, cols, d, after):
    ws = wb.create_sheet(f"{data_title} - columns", wb.sheetnames.index(after) + 1)
    ws.append(["column", "source", "what it is"])
    for c in cols:
        src, txt = d.get(c, ("", "(undocumented)"))
        ws.append([c, src, txt])
    for c in ws[1]: c.font = Font(bold=True)
    ws.freeze_panes = "A2"
    ws.column_dimensions["A"].width = 30; ws.column_dimensions["B"].width = 38; ws.column_dimensions["C"].width = 120
    for row in ws.iter_rows(min_row=2):
        for c in row: c.alignment = Alignment(wrap_text=True, vertical="top")
    return [[ws.cell(r, k).value for k in (1, 2, 3)] for r in range(1, ws.max_row + 1)]

wb = Workbook(); about = wb.active; about.title = "About"
lines = [
    ("Cyp / TE pipeline results: Atallah lab pipeline (archival reconstruction), full run 2026-11-22/23", True),
    ("Only this pipeline's outputs are in this workbook. Lab scripts were run unchanged (arguments and flags only).", False),
    ("", False),
    ("Tabs", True),
    ("Per species: one row per species that ran stage 3 (301). Columns by stage:", False),
    ("  ecology: exposure_group (high/low/none) from the ecology sheet's 'Insecticide exposure' column; 'none' = Possible/Unknown, left out of the tests.", False),
    ("  confounder: low species whose detox genes track host toxins (cactophiles, wild-fungus breeders, D. sechellia).", False),
    ("  stage 1 (RepeatModeler/RepeatMasker .tbl): genome_bp, pct_masked, pct_LINE/LTR/DNA/RC/unclassified/simple, repeat_families.", False),
    ("  stage 3 (repeatOpp/Locate_TE/CleanAnnasse tables, window gene +/-3 kb): cyp_loci, te_hits (TE classes), simple_lowcomplex_hits, te_per_cyp_locus, xeno_*.", False),
    ("  s4_<test>_*: stage 4 per-species summary from each compare script (blank for 'none' species, which are not in the config).", False),
    ("Each data tab has a '<tab> - columns' tab after it describing every column.", False),
    ("Statistics: the three compare scripts' tests, copied from their reports.", False),
    ("Cyp genes (stage 4): every Cyp gene the exposure script counted (high/low species only).", False),
    ("Xenobiotic genes (stage 4): the 10-gene list's matches per species.", False),
    ("", False),
    ("Read before interpreting", True),
    ("1. Every Cyp gene has >=1 'TE', so Fisher/chi-square have no variance (p = 1 is saturation, not a null result). The stage 3 tables only contain genes with a hit, "
     "and the counts include Simple_repeat / Low_complexity hits, not just TEs.", False),
    ("2. Mann-Whitney, Welch and the species bootstrap all go the OPPOSITE way from the hypothesis: low-exposure species have ~1.4x more hits per Cyp gene.", False),
    ("3. Report text says '5 species split 3-2' in places: that is fixed wording in the lab scripts, not this run's numbers (235 species).", False),
    ("4. Gene-level tests are pseudoreplicated (genes within a species are not independent). The species bootstrap is the more honest estimate.", False),
    ("5. Known lab-script defects are preserved: CleanAnnasse substring name matching (e.g. Cyp12d1 named but absent from the annotation), repeatOpp duplicate rows.", False),
    ("6. Exposure groups and the xenobiotic gene list are not lab records; they were built for this run (ecology sheet rule; list citations still to be checked).", False),
    ("7. Median ~21 Cyp loci per species vs ~85 real Cyp genes in Drosophila: renaming/annotation coverage limits what is counted.", False),
    ("", False),
    ("Files: /bioinformatics/atallah-repeats (stage folders 1-4, summary/). Stage 4 reports and plots: 4-combine-and-compare/output/compare/.", False),
]
for t, b in lines:
    about.append([t]); about.cell(about.max_row, 1).font = Font(bold=b)
about.column_dimensions["A"].width = 160

sheet(wb, "Per species", sp)
ws = sheet(wb, "Statistics", st)
for row in ws.iter_rows(min_row=2):
    for c in row: c.alignment = Alignment(wrap_text=True, vertical="top")
sheet(wb, "Cyp genes (stage 4)", tsv(f"{C}/exposure.per-gene.csv", ","))
sheet(wb, "Xenobiotic genes (stage 4)", tsv(f"{C}/xenobiotic.per-gene.csv", ","))
import json
dump = {}
for title, d in [("Per species", PER_SPECIES), ("Statistics", STATS), ("Cyp genes (stage 4)", PER_GENE), ("Xenobiotic genes (stage 4)", PER_GENE)]:
    cols = [c.value for c in wb[title][1]]
    dump[f"{title} - columns"] = doc_tab(wb, title, cols, d, title)
    assert all(r[2] != "(undocumented)" for r in dump[f"{title} - columns"]), title
json.dump(dump, open(f"{A}/summary/column-docs.json", "w"), indent=1)
out = f"{A}/summary/pipeline-results.xlsx"; wb.save(out)
for s in st: print(s)
print("wrote", out)
