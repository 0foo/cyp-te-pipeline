# NEW (not a lab file, 2026-11-23): builds summary/species-hits-and-stage4.xlsx (the Google Sheet
# "Cyp-TE per-species hits and stage 4"), a workbook with one row per species (301).
# COLUMNS: the stage 3 totals from summary/per-species.tsv (Cyp loci, repeat hits within gene +/- 3 kb,
#   how many are TEs, simple and low-complexity hits, TEs per locus). Then, for each stage 4 compare
#   script (all / cncc / xeno), its per-species summary from output/compare/<s>.species-summary.csv:
#   n_cyp_genes, genes with a TE, total "TEs", TEs per kb and TEs per gene, plus the TFBS count for cncc.
# Nothing is recomputed from the stage 4 outputs. The compare scripts' tests are pooled across species
#   (there is no per-species test), so the per-species values shown are the inputs to those tests.
#   Stage 4 columns are blank for the 66 "none" exposure species.
# TABS: About, Per species, Columns (a definition of every column).
# Usage: nix-shell -p 'python3.withPackages(p: [p.openpyxl])' --run 'python3 -B make-hits-sheet.py'
import csv
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment
from openpyxl.utils import get_column_letter

A = "/bioinformatics/atallah-repeats"
C = f"{A}/4-combine-and-compare/output/compare"
tsv = lambda p, d="\t": list(csv.DictReader(open(p), delimiter=d))
sp = tsv(f"{A}/summary/per-species.tsv")
s4 = {t: {r["species"]: r for r in tsv(f"{C}/{t}.species-summary.csv", ",")} for t in ("exposure", "cncc", "xenobiotic")}

def num(v):
    if v in ("", None): return ""
    try:
        f = float(v); return int(f) if f.is_integer() and "." not in str(v) else round(f, 3)
    except ValueError:
        return v

COLS = [  # (column, source, description, getter)
    ("species", "ID", "Species label used in this pipeline.", lambda r, e, c, x: r["species"]),
    ("exposure_group", "Ecology sheet", "high / low / none. 'none' (Possible, Unknown) species are not in the stage 4 config, so their stage 4 columns are blank.", lambda r, e, c, x: r["exposure_group"]),
    ("confounder", "make-exposure.py", "'host-toxin confounder' for low species whose detox genes track plant/fungal toxins; blank otherwise.", lambda r, e, c, x: r["confounder"]),
    ("in_stage4_tests", "species_config.ini", "yes if the species was in the stage 4 comparisons (high or low exposure).", lambda r, e, c, x: "yes" if e else "no"),
    ("cyp_loci", "Stage 3", "Distinct Cyp loci repeatOpp.py kept.", lambda r, e, c, x: r["cyp_loci"]),
    ("total_hits", "Stage 3", "Total RepeatMasker hits in the Cyp windows (gene +/-3 kb, hit fully inside), de-duplicated: TEs plus simple/low-complexity repeats.", lambda r, e, c, x: r["unique_hits"]),
    ("total_TEs", "Stage 3", "Of total_hits, those that are transposable elements (LINE, LTR, DNA, RC, SINE, Unknown...; excludes Simple_repeat, Low_complexity, Satellite, RNA, ARTEFACT).", lambda r, e, c, x: r["te_hits"]),
    ("simple_lowcomplex_hits", "Stage 3", "total_hits - total_TEs.", lambda r, e, c, x: r["simple_lowcomplex_hits"]),
    ("pct_hits_TE", "Stage 3", "total_TEs / total_hits x 100.", lambda r, e, c, x: r["pct_hits_te"]),
    ("TEs_per_cyp_locus", "Stage 3", "total_TEs / cyp_loci.", lambda r, e, c, x: r["te_per_cyp_locus"]),
]
D = {"n_cyp_genes": "Cyp genes the script counted (unique Name= in the combined GFF3)",
     "n_genes_with_te": "of those, genes with >=1 'TE'",
     "pct_genes_with_te": "% of counted genes with >=1 'TE'",
     "total_te_count": "total 'TEs' across those genes. NOTE: the scripts count every repeat hit as a TE (simple repeats included)",
     "te_per_kb": "total_te_count per kb of Cyp gene length",
     "te_per_gene": "total_te_count / n_cyp_genes: the per-species value behind the pooled tests"}
for t, lab, what in [("exposure", "all", "all Cyp genes"), ("cncc", "cncc", "Cyp genes near a predicted CncC:Maf-S site"),
                     ("xenobiotic", "xeno", "the 10 xenobiotic-list genes")]:
    if t == "cncc":
        COLS.append(("cncc_tf_binding_sites", "Stage 4: compare_te_cyp_cncc.py", "All TF binding sites FIMO found in the species (0 = motif scan failed).",
                     lambda r, e, c, x: x["cncc"].get("n_tfbs_total", "")))
    for k, d in D.items():
        COLS.append((f"{lab}_{k}", f"Stage 4: compare_te_cyp_{t}.py", f"{what}: {d}.",
                     (lambda t, k: lambda r, e, c, x: x[t].get(k, ""))(t, k)))

wb = Workbook(); ab = wb.active; ab.title = "About"
for line, bold in [
    ("Per-species hits, TEs and stage 4 results: Cyp / TE pipeline, run 2026-11-22/23", True),
    ("One row per species (301). Only this pipeline's outputs. Lab scripts run unchanged.", False),
    ("", False),
    ("total_hits / total_TEs: from stage 3 (RepeatMasker hits within +/-3 kb of each Cyp gene).", False),
    ("all_* / cncc_* / xeno_*: each stage 4 compare script's per-species output, for its gene set.", False),
    ("The stage 4 tests (Fisher, Mann-Whitney, Welch, species bootstrap) are pooled across species, so there is no per-species p-value; "
     "the per-species columns here are what the tests are built from. Pooled results: all three 'Contradict the hypothesis' (see the main results sheet).", False),
    ("Stage 4 'TE' counts include simple/low-complexity repeats (every stage 3 hit), so compare them with total_hits, not total_TEs. "
     "They can be a little lower than total_hits because the scripts key genes by name (hits on loci whose name is missing from the annotation drop out).", False),
    ("Stage 4 columns are blank for the 66 'none' exposure species (not in the test config).", False),
    ("The 'Columns' tab explains every column.", False),
]:
    ab.append([line]); ab.cell(ab.max_row, 1).font = Font(bold=bold)
ab.column_dimensions["A"].width = 160

ws = wb.create_sheet("Per species")
ws.append([c[0] for c in COLS])
for r in sp:
    e = r["species"] in s4["exposure"]
    ws.append([num(g(r, e, None, {t: s4[t].get(r["species"], {}) for t in s4})) for _, _, _, g in COLS])
for c in ws[1]: c.font = Font(bold=True)
ws.freeze_panes = "B2"; ws.auto_filter.ref = ws.dimensions
for i, c in enumerate(COLS, 1): ws.column_dimensions[get_column_letter(i)].width = max(10, len(c[0]) + 2)

cw = wb.create_sheet("Columns")
cw.append(["column", "source", "what it is"])
for n, s, d, _ in COLS: cw.append([n, s, d])
for c in cw[1]: c.font = Font(bold=True)
cw.freeze_panes = "A2"
for col, w in zip("ABC", (28, 36, 130)): cw.column_dimensions[col].width = w
for row in cw.iter_rows(min_row=2):
    for c in row: c.alignment = Alignment(wrap_text=True, vertical="top")
out = f"{A}/summary/species-hits-and-stage4.xlsx"; wb.save(out)
print("wrote", out, len(sp), "rows,", len(COLS), "columns")
