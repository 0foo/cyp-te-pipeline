# NEW (not a lab file, 2026-11-24): builds cyp-te-validation.xlsx (one row per species, 301) for upload as a
# Google Sheet. Counts only what the pipeline's scripts produced:
#   Cyp tab: the Cyp loci repeatOpp.py selected (runs/<LABEL>/filtered.gff) and the names CleanAnnasse.py
#            gave them (runs/<LABEL>/DAnasse_TE_Cyp.txt column 2). Each locus is then checked: is it really a
#            cytochrome P450, or something sharing the 'Cyp' prefix (cyclophilin Cyp1/Cyp40/Cypl/Cyp33,
#            P450 reductase Cpr, ...)?
#   Transposons tab: the species' RepeatMasker output (stage 1, rmodeler/out/<rm_sample>.rm.out), genome-wide,
#            plus the TEs stage 3 found in the Cyp windows (summary/per-species.tsv).
# Usage: python3 -B make-validation-sheet.py   (needs openpyxl)
import csv, re, collections
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment
from openpyxl.utils import get_column_letter
import sys; sys.path.insert(0, __import__('os').path.dirname(__file__)); from column_labels import LABELS, TE_KEEP

A = "/bioinformatics/atallah-repeats"
R = f"{A}/3-pair-te-with-cyp/runs"
P450 = re.compile(r"^Cyp\d+[a-z]+\d*", re.I)   # P450 names: Cyp4g1, Cyp12d1-p, Cyp313a, CYP9F2 (not Cyp1/Cyp40/Cypl)
HALLOWEEN = {"phm", "dib", "sad", "shd", "spo", "spok"}      # P450s FlyBase names by symbol
NONTE = {"Simple_repeat", "Low_complexity", "Satellite", "rRNA", "tRNA", "snRNA", "scRNA", "srpRNA", "ARTEFACT"}
TOPS = ["LINE", "LTR", "DNA", "RC", "SINE", "Unknown"]
REF = 84   # P450 genes per Drosophila genome: Good et al. 2014 found 74-94 across 12 species; 84 = middle

HALL_NUM = {"Cyp306a1": "phm", "Cyp302a1": "dib", "Cyp315a1": "sad", "Cyp314a1": "shd", "Cyp307a1": "spo", "Cyp307a2": "spok"}
# Scored: the 4 Halloween genes Good et al. 2014 (GBE 6:1118, Table 1) found single-copy in all 12 species.
# spook/spookier (Cyp307a1/a2) are 'unstable' there (a loss / a duplication), so they are not expected in every species.
HALL_STABLE = ["Cyp306a1", "Cyp302a1", "Cyp315a1", "Cyp314a1"]
HALL_SYM = {v: k for k, v in HALL_NUM.items()}
def halloween(path, genes_only):
    """Halloween P450s (by number or FlyBase symbol) named in a GFF's Name= values."""
    found = set()
    for l in open(path):
        x = l.rstrip("\n").split("\t")
        if len(x) < 9 or (genes_only and x[2] != "gene"): continue
        for n in attr(x[8], "Name").split(","):
            if n in HALL_NUM: found.add(n)
            elif n in HALL_SYM: found.add(HALL_SYM[n])
    return found
def gff_for(lab):
    return (f"{A}/2-rename-genes/input/gffs/DROSOPHILA_MELANOGASTER_final.gff" if lab == "DROSOPHILA_MELANOGASTER"
            else f"{A}/3-pair-te-with-cyp/input/renamed-gff/{lab}_final_withDmelNames.gff")
# Conserved single-copy Cyps with a known conserved function: 'stable' clades of Good et al. 2014 (Table 1; one gene in
# each of 12 Drosophila genomes) that the paper assigns a function (Results, "P450 Gene Gain and Loss among Species"):
# ecdysteroid synthesis Cyp302a1/dib, Cyp314a1/shd, Cyp306a1/phm; ecdysone modification Cyp18a1; bristle development
# Cyp303a1; cuticular hydrocarbon metabolism Cyp4g1. Plus Cyp315a1/sad: stable in Table 1, ecdysone C2-hydroxylase
# (Warren et al. 2002, PNAS 99:11043), not named in that sentence.
CONSERVED = {"Cyp18a1": "Cyp18a1", "Cyp303a1": "Cyp303a1", "Cyp4g1": "Cyp4g1",
             "Cyp302a1": "dib", "Cyp314a1": "shd", "Cyp306a1": "phm", "Cyp315a1": "sad"}
CONS_ALIAS = {**{k: k for k in CONSERVED}, **{v: k for k, v in CONSERVED.items()}}
def conserved(path, genes_only):
    found = set()
    for l in open(path):
        x = l.rstrip("\n").split("\t")
        if len(x) < 9 or (genes_only and x[2] != "gene"): continue
        found.update(CONS_ALIAS[n] for n in attr(x[8], "Name").split(",") if n in CONS_ALIAS)
    return found
def is_p450(name): return bool(P450.match(name)) or name in HALLOWEEN
def attr(a, k):
    m = re.search(rf"(?:^|;){k}=([^;]*)", a); return m.group(1) if m else ""

sp = list(csv.DictReader(open(f"{A}/summary/per-species.tsv"), delimiter="\t"))

# ---- Cyp
cyp_rows = []
for r in sp:
    lab = r["species"]; loci = {}
    for l in open(f"{R}/{lab}/filtered.gff"):
        x = l.rstrip("\n").split("\t")
        if len(x) >= 9: loci.setdefault(attr(x[8], "ID") or x[8], [n for n in attr(x[8], "Name").split(",") if n])
    final = {x[1] for x in (l.rstrip("\n").split("\t") for l in open(f"{R}/{lab}/DAnasse_TE_Cyp.txt")) if len(x) > 2}
    non = [",".join(n) or "(no Name)" for n in loci.values() if not any(is_p450(k) for k in n)]
    non_final = sorted(n for n in final if not is_p450(n))
    cyp_rows.append(dict(species=lab, exposure_group=r["exposure_group"], total_cyp_count=len(loci), unique_cyp_count=len(final),
                         not_p450_loci=len(non), pct_not_p450=round(100 * len(non) / len(loci), 1) if loci else 0,
                         not_p450_names="; ".join(sorted(set(non) | set(non_final))),
                         pct_absent=round(max(0, 100 * (1 - len(final) / REF)), 1)))
    hf = halloween(f"{R}/{lab}/filtered.gff", False) & set(HALL_STABLE)
    ha = halloween(gff_for(lab), True) & set(HALL_STABLE)
    cyp_rows[-1].update(percentage_halloween_found=round(100 * len(hf) / len(HALL_STABLE), 1),
                        halloween_found=", ".join(f"{HALL_NUM[n]} ({n})" for n in HALL_STABLE if n in hf),
                        percentage_halloween_in_annotation=round(100 * len(ha) / len(HALL_STABLE), 1))
    cf, ca = conserved(f"{R}/{lab}/filtered.gff", False), conserved(gff_for(lab), True)
    cyp_rows[-1].update(percentage_conserved_cyps_found=round(100 * len(cf) / len(CONSERVED), 1),
                        conserved_cyps_missing=", ".join(k if CONSERVED[k] == k else f"{CONSERVED[k]} ({k})" for k in CONSERVED if k not in cf),
                        percentage_conserved_cyps_in_annotation=round(100 * len(ca) / len(CONSERVED), 1))

# ---- Transposons (genome-wide, RepeatMasker .out)
loci_by_sp = collections.defaultdict(list)
for l in csv.DictReader(open(f"{A}/summary/cyp-genes.tsv"), delimiter="\t"): loci_by_sp[l["species"]].append(l)
te_rows = []
for r in sp:
    rows = te = te_bp = nonte = simple = 0; ids = set(); top = collections.Counter(); fam = set()
    with open(f"{A}/rmodeler/out/{r['rm_sample']}.rm.out") as fh:
        for i, l in enumerate(fh):
            x = l.split()
            if i < 3 or len(x) < 15: continue
            rows += 1
            cls = x[10]; t = cls.split("/")[0]
            if t in NONTE:
                nonte += 1; simple += t == "Simple_repeat"; continue
            te += 1; te_bp += int(x[6]) - int(x[5]) + 1
            ids.add(x[14]); fam.add(x[9])
            top[t if t in TOPS else "other"] += 1
    gb = int(r["genome_bp"])
    row = dict(species=r["species"], exposure_group=r["exposure_group"], rm_sample=r["rm_sample"],
               total_te_count=len(ids), te_fragments=te, repeatmasker_rows=rows, simple_lowcomplex_rows=nonte,
               pct_not_true_te=round(100 * nonte / rows, 1) if rows else 0,
               pct_microsatellite=round(100 * simple / rows, 1) if rows else 0,
               pct_not_true_te_near_cyp=round(100 - float(r["pct_hits_te"]), 1),
               te_bp=te_bp, pct_genome_te=round(100 * te_bp / gb, 2), te_per_mb=round(len(ids) / (gb / 1e6), 1),
               te_families_hit=len(fam), repeat_families_in_library=int(r["repeat_families"]), pct_masked_rm_tbl=float(r["pct_masked"]))
    for t in TOPS + ["other"]: row[f"te_{t}"] = top[t]
    row.update(te_near_cyp=int(r["te_hits"]), te_near_cyp_per_locus=float(r["te_per_cyp_locus"]))
    # stage 4: the 'TEs' compare_te_cyp_exposure.py counts = mobile_genetic_element rows whose
    # Description "Within range of <gene>" names a gene feature in the combined GFF3; split by repeat_class
    genes, s4 = set(), []
    for l in open(f"{A}/4-combine-and-compare/output/combined-gff/{r['species']}.combined.gff3"):
        if l.startswith("#"): continue
        x = l.rstrip("\n").split("\t")
        if len(x) < 9: continue
        if x[2] == "gene": genes.add(attr(x[8], "Name"))
        m = re.search(r"Within range of (\S+)", attr(x[8], "Description"))
        if m: s4.append((m.group(1), attr(x[8], "repeat_class").split("/")[0]))
    s4 = [c for g, c in s4 if g in genes]
    s4_non = sum(c in NONTE for c in s4)
    row.update(s4_te_count=len(s4), s4_true_te=len(s4) - s4_non,
               s4_pct_not_true_te=round(100 * s4_non / len(s4), 1) if s4 else "",
               s4_pct_microsatellite=round(100 * sum(c == "Simple_repeat" for c in s4) / len(s4), 1) if s4 else "")
    kb = sum(int(l["end"]) - int(l["start"]) + 1 + 6000 for l in loci_by_sp[r["species"]]) / 1000
    near = int(r["te_hits"]) / kb if kb else 0
    gw = te / (gb / 1000)
    row.update(te_per_10kb_near_cyp=round(10 * near, 3), te_per_10kb_genome=round(10 * gw, 3),
               near_cyp_vs_genome_te_density=round(near / gw, 3) if gw else "")
    te_rows.append(row)
    print(r["species"], row["total_te_count"], flush=True)

CYP_DOC = {
    "species": "Species label used in this pipeline.",
    "exposure_group": "high / low / none (ecology sheet; none = not in the stage 4 tests).",
    "total_cyp_count": "Cyp loci repeatOpp.py selected from the renamed annotation (distinct gene IDs in filtered.gff; repeatOpp's duplicate rows counted once).",
    "unique_cyp_count": "Distinct Cyp names after CleanAnnasse.py (column 2 of DAnasse_TE_Cyp.txt). Lower than total_cyp_count because paralogs and multi-name loci collapse to one name, including the substring-matching defect (Cyp4g15 -> Cyp4g1).",
    "not_p450_loci": "Loci the scripts selected that are not cytochrome P450s: no name of the form Cyp<number><letters><number> or a Halloween symbol (phm, dib, sad, shd, spo, spok), e.g. cyclophilins (Cyp1, Cyp40, Cypl, Cyp33) or P450 reductase (Cpr).",
    "pct_not_p450": "not_p450_loci / total_cyp_count x 100.",
    "not_p450_names": "Names of any non-P450 loci or final names (blank = none).",
    "percentage_halloween_found": "Share of the 4 stable Halloween P450s (phm/Cyp306a1, dib/Cyp302a1, sad/Cyp315a1, shd/Cyp314a1; ecdysone synthesis) among the loci the scripts selected (repeatOpp.py filtered.gff), as %. These 4 are single-copy in all 12 Drosophila genomes in Good et al. 2014 (GBE 6:1118, Table 1); spook and spookier are left out because that paper finds a loss and a duplication in them. 0 in every species except D. melanogaster: the lab gene list names them by number (Cyp306a1...) but the renamed annotations name them by symbol (phm...), so the scripts' match misses them. D. melanogaster is found only because its RefSeq records list the numbers as synonyms.",
    "halloween_found": "Which of the 4 the scripts selected.",
    "percentage_halloween_in_annotation": "Share of the 4 present in the species' renamed annotation (stage 2 output, the file the scripts search), as %. The gap to percentage_halloween_found is genes lost at selection; below 100% here means lost at annotation/renaming.",
    "percentage_conserved_cyps_found": "Share of 7 conserved single-copy Cyps with known conserved functions among the loci the scripts selected, as %. The 7 are 'stable' (one copy in each of 12 Drosophila genomes) in Good et al. 2014, Table 1, with functions from that paper: ecdysteroid synthesis dib/Cyp302a1, shd/Cyp314a1, phm/Cyp306a1; ecdysone modification Cyp18a1; bristle development Cyp303a1; cuticular hydrocarbon metabolism Cyp4g1; plus sad/Cyp315a1 (ecdysone synthesis; Warren et al. 2002). A missing one points to acquisition failure. The scripts can reach at most 3/7 (42.9%) outside D. melanogaster: the 4 Halloween genes are named by symbol in the annotations.",
    "conserved_cyps_missing": "Which of the 7 the scripts did not select.",
    "percentage_conserved_cyps_in_annotation": "Share of the 7 present in the species' renamed annotation (the file the scripts search), as %. Below 100% = lost before the scripts ran (annotation or renaming).",
    "pct_absent": f"1 - unique_cyp_count / {REF}, as %, floored at 0. {REF} = middle of the 74-94 P450 genes per genome Good et al. 2014 found across 12 Drosophila species.",
}
TE_DOC = {
    "species": "Species label used in this pipeline.",
    "exposure_group": "high / low / none.",
    "rm_sample": "RepeatMasker sample (which assembly).",
    "total_te_count": "Transposable elements RepeatMasker identified genome-wide: distinct element IDs (.out column 15), so fragments of one insertion count once. Excludes Simple_repeat, Low_complexity, Satellite, small RNAs, ARTEFACT.",
    "te_fragments": "TE rows in the .out file (fragments counted separately).",
    "repeatmasker_rows": "All rows in the .out file (TEs + simple/low-complexity + other).",
    "simple_lowcomplex_rows": "Non-TE rows (simple, low-complexity, satellite, RNA, artefact).",
    "pct_not_true_te": "simple_lowcomplex_rows / repeatmasker_rows x 100: share of everything RepeatMasker reported genome-wide that is not a transposon (microsatellites, low-complexity, satellites, small RNAs, artefacts).",
    "pct_microsatellite": "Simple_repeat rows (microsatellites, e.g. (TGTC)n) / repeatmasker_rows x 100. The part of pct_not_true_te that is microsatellites alone.",
    "pct_not_true_te_near_cyp": "Same as pct_not_true_te but for the stage 3 hits in the Cyp windows (100 - pct of hits that are TEs). These non-TE hits are still counted as 'TEs' by the stage 4 compare scripts.",
    "te_bp": "Sum of TE row lengths, bp. Overlapping rows are counted twice, so this can slightly exceed the masked bp (pct_genome_te can be above pct_masked_rm_tbl in TE-rich genomes).",
    "pct_genome_te": "te_bp / genome length x 100.",
    "te_per_mb": "total_te_count per Mb of assembly.",
    "te_families_hit": "Distinct repeat families with >=1 TE hit.",
    "repeat_families_in_library": "Families in the species' RepeatModeler library (families.fa).",
    "pct_masked_rm_tbl": "% of genome masked, all repeat types, from RepeatMasker's own .tbl summary (cross-check for pct_genome_te).",
    **{f"te_{t}": f"TE rows of class {t}." for t in TOPS}, "te_other": "TE rows of any other class (e.g. Retroposon, PLE).",
    "te_near_cyp": "TEs stage 3 found within +/-3 kb of the Cyp loci (de-duplicated; simple repeats excluded).",
    "te_near_cyp_per_locus": "te_near_cyp / Cyp loci.",
    "s4_te_count": "Hits the stage 4 compare scripts count as 'TEs' for this species' Cyp genes (combined GFF3 mobile_genetic_element rows 'Within range of' a Cyp gene; same rule as compare_te_cyp_exposure.py, so it equals all_total_te_count in the results sheet for high/low species). Filled for all 301 species, since the combined GFF3 was built for all.",
    "s4_true_te": "Of s4_te_count, rows whose repeat_class is a transposon.",
    "s4_pct_not_true_te": "Share of s4_te_count that is not a transposon (Simple_repeat, Low_complexity, Satellite, small RNAs, ARTEFACT), as %: how much of what stage 4 tested is microsatellite/low-complexity sequence.",
    "s4_pct_microsatellite": "Share of s4_te_count that is Simple_repeat (microsatellites), as %.",
    "te_per_10kb_near_cyp": "te_near_cyp per 10 kb of Cyp window (window = locus length + 3 kb each side; overlapping windows of clustered Cyps count once per locus, as stage 3 counts their TEs).",
    "te_per_10kb_genome": "te_fragments per 10 kb of assembly (genome-wide; same unit as stage 3, one RepeatMasker row per hit).",
    "near_cyp_vs_genome_te_density": "te_per_10kb_near_cyp / te_per_10kb_genome. 1 = Cyp regions carry TEs at the genome average; >1 enriched, <1 depleted. Controls for each genome's overall TE load. Biased low: stage 3 keeps only TEs lying fully inside the window, and gene-rich euchromatin is TE-poor anyway.",
}

wb = Workbook(); ab = wb.active; ab.title = "About"
for t, b in [("Cyp and transposon validation: Cyp / TE pipeline (2026-11-24)", True),
             ("One row per species (301). Only what this pipeline's scripts produced; nothing re-annotated.", False),
             ("Cyp: the loci the lab scripts selected and named (repeatOpp.py, CleanAnnasse.py), checked for non-P450 genes.", False),
             ("Transposons: RepeatMasker output genome-wide (stage 1) and near the Cyp loci (stage 3).", False),
             ("Each data tab has a '- columns' tab after it explaining every column.", False)]:
    ab.append([t]); ab.cell(ab.max_row, 1).font = Font(bold=b)
ab.column_dimensions["A"].width = 120

def add(title, rows, doc, cols=None):
    ws = wb.create_sheet(title); cols = cols or list(rows[0])
    ws.append([LABELS.get(c, c) for c in cols])
    for r in rows: ws.append([r[c] for c in cols])
    for c in ws[1]: c.font = Font(bold=True); c.alignment = Alignment(wrap_text=True, vertical="top")
    ws.freeze_panes = "B2"; ws.auto_filter.ref = ws.dimensions
    for i, c in enumerate(cols, 1): ws.column_dimensions[get_column_letter(i)].width = 22
    ws.column_dimensions["A"].width = 34
    cw = wb.create_sheet(f"{title} - columns"); cw.append(["column", "what it is", "short name"])
    for c in cols: cw.append([LABELS.get(c, c), doc[c], c])
    for c in cw[1]: c.font = Font(bold=True)
    cw.column_dimensions["A"].width = 28; cw.column_dimensions["B"].width = 140
    for row in cw.iter_rows(min_row=2):
        for c in row: c.alignment = Alignment(wrap_text=True, vertical="top")
add("Cyp", cyp_rows, CYP_DOC)
add("Transposons", te_rows, TE_DOC, TE_KEEP)
out = f"{A}/summary/cyp-te-validation.xlsx"; wb.save(out); print("wrote", out)
