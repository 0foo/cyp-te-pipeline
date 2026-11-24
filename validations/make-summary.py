# NEW (not a lab file, 2026-11-23): per-species and per-Cyp-locus summary of the cyp-te-pipeline run
# in /bioinformatics/atallah-repeats. It reads only that run's files and writes three TSVs to summary/.
# The two workbook builders, make-validation-sheet.py and sheets/make-*.py, read these TSVs.
#
# INPUTS
#   3-pair-te-with-cyp/pairs.tsv (rows with status RUN: 301 species) and runs/<LABEL>/filtered.gff
#   and GenesAffectedByTEs.txt (stage 3); rmodeler/out/<rm_sample>.rm.tbl and -families.fa (stage 1);
#   4-combine-and-compare/exposure.tsv and xenobiotic_resistance_cyp_genes.txt (stage 4 settings).
#
# HOW IT WORKS (per species)
#   - Cyp loci: the filtered.gff lines de-duplicated by gene ID, which undoes repeatOpp.py's duplicate rows.
#   - Hits: GenesAffectedByTEs.txt rows de-duplicated on (locus ID, RepeatMasker row). That file
#     is used, rather than DAnasse_TE_Cyp.txt, because it still carries the locus ID. Each hit is
#     classed by the top level of its class/family: "nonTE" for Simple_repeat, Low_complexity,
#     Satellite, rRNA, tRNA, snRNA, scRNA, srpRNA and ARTEFACT, and "TE" otherwise (LINE, LTR, DNA,
#     RC, SINE, Unknown, ...).
#   - Window size: gene length + 6,000 bp per locus (+/- 3 kb). That gives TE hits per 10 kb of window.
#   - Genome: size, % masked and per-class % from the .rm.tbl; library size = number of ">" lines in -families.fa.
#   - Xenobiotic targets: a target is present when it equals (ignoring case) one of a locus's Name= values.
# OUTPUTS (summary/, tab-separated)
#   per-species.tsv       one row per species (301): exposure fields, genome stats, Cyp loci, hits, TE/nonTE
#   xenobiotic-genes.tsv  one row per species x target gene: present, loci, TE and non-TE hits
#   cyp-genes.tsv         one row per species x Cyp locus: names, coordinates, hits by class
# Usage: python3 -B make-summary.py   (standard library only; writes into /bioinformatics/atallah-repeats/summary)
import csv, os, re, collections
A = "/bioinformatics/atallah-repeats"
OUT = f"{A}/summary"
NONTE = {"Simple_repeat", "Low_complexity", "Satellite", "rRNA", "tRNA", "snRNA", "scRNA", "srpRNA", "ARTEFACT"}
TOPS = ["LINE", "LTR", "DNA", "RC", "SINE", "Unknown"]

pairs = {r["gff_label"]: r for r in csv.DictReader(open(f"{A}/3-pair-te-with-cyp/pairs.tsv"), delimiter="\t") if r["status"] == "RUN"}
expo = {r["gff_label"]: r for r in csv.DictReader(open(f"{A}/4-combine-and-compare/exposure.tsv"), delimiter="\t")}
targets = [l.split("#")[0].strip() for l in open(f"{A}/4-combine-and-compare/xenobiotic_resistance_cyp_genes.txt")]
targets = [t for t in targets if t]

def tbl(s):
    t = open(f"{A}/rmodeler/out/{s}.rm.tbl").read()
    g = lambda pat: (float(m.group(1)) if (m := re.search(pat, t)) else None)
    return dict(genome_bp=int(g(r"total length:\s+(\d+) bp")), genome_seqs=int(g(r"sequences:\s+(\d+)")),
                pct_masked=g(r"bases masked:.*?\(\s*([\d.]+) %\)"), pct_LINE=g(r"LINEs:\s+\d+\s+\d+ bp\s+([\d.]+) %"),
                pct_LTR=g(r"LTR elements:\s+\d+\s+\d+ bp\s+([\d.]+) %"), pct_DNA=g(r"DNA transposons\s+\d+\s+\d+ bp\s+([\d.]+) %"),
                pct_RC=g(r"Rolling-circles\s+\d+\s+\d+ bp\s+([\d.]+) %"), pct_unclassified=g(r"Unclassified:\s+\d+\s+\d+ bp\s+([\d.]+) %"),
                pct_simple=g(r"Simple repeats:\s+\d+\s+\d+ bp\s+([\d.]+) %"))

def names(attr):
    m = re.search(r"(?:^|;)Name=([^;]*)", attr)
    return m.group(1).split(",") if m else []

def gid(attr):
    m = re.search(r"(?:^|;)ID=([^;]*)", attr)
    return m.group(1) if m else attr

sp_rows, xg_rows, cg_rows = [], [], []
for lab in sorted(pairs):
    p, e = pairs[lab], expo.get(lab, {})
    s = p["rm_sample"]; run = f"{A}/3-pair-te-with-cyp/runs/{lab}"
    nfam = sum(1 for l in open(f"{A}/rmodeler/out/{s}-families.fa") if l.startswith(">"))
    loci = {}                      # Cyp loci kept by repeatOpp (filtered.gff), deduplicated by ID
    raw_loci = 0
    for l in open(f"{run}/filtered.gff"):
        x = l.rstrip("\n").split("\t")
        if len(x) < 9: continue
        raw_loci += 1
        loci.setdefault(gid(x[8]), dict(seq=x[0], start=int(x[3]), end=int(x[4]), names=names(x[8])))
    rows = 0; hits = {}            # unique (locus, RepeatMasker row)
    # GenesAffectedByTEs.txt = DAnasse_TE_Cyp.txt line for line, but keeps the GFF attributes (locus ID);
    # CleanAnnasse.py only replaces column 2 with a gene name.
    for l in open(f"{run}/GenesAffectedByTEs.txt"):
        x = l.rstrip("\n").split("\t")
        if len(x) < 3: continue
        rows += 1
        f = x[2].split()
        hits.setdefault((gid(x[1]), x[2].strip()), dict(cls=f[10], rep=f[9]))
    per_locus = collections.defaultdict(lambda: collections.Counter())
    for (g, _), h in hits.items():
        top = h["cls"].split("/")[0]
        kind = "nonTE" if top in NONTE else "TE"
        per_locus[g][kind] += 1
        per_locus[g]["TE_" + top if kind == "TE" else "nonTE_" + top] += 1
    te_total = sum(c["TE"] for c in per_locus.values())
    nonte_total = sum(c["nonTE"] for c in per_locus.values())
    loci_te = sum(1 for g in loci if per_locus[g]["TE"] > 0)
    loci_any = sum(1 for g in loci if per_locus[g]["TE"] + per_locus[g]["nonTE"] > 0)
    by_top = collections.Counter()
    for c in per_locus.values():
        for k, v in c.items():
            if k.startswith("TE_"): by_top[k[3:]] += v
    kb = sum(v["end"] - v["start"] + 1 + 6000 for v in loci.values()) / 1000
    row = dict(species=lab, rm_sample=s, sheet_label=e.get("sheet_label", ""), exposure_group=e.get("exposure", ""),
               insecticide_exposure=e.get("insecticide_exposure", ""), exposure_setting=e.get("exposure_setting", ""),
               resistance_evidence=e.get("resistance_evidence", ""), confounder=e.get("confounder", ""))
    row.update(tbl(s)); row["repeat_families"] = nfam
    row.update(gene_frac_on_rm_seqs=p["gene_frac_on_rm_contigs"], cyp_genes_in_annotation=p["cyp_genes"],
               cyp_loci=len(loci), cyp_loci_rows_raw=raw_loci, te_table_rows_raw=rows, unique_hits=len(hits),
               te_hits=te_total, simple_lowcomplex_hits=nonte_total,
               pct_hits_te=round(100 * te_total / len(hits), 1) if hits else 0,
               cyp_loci_with_te=loci_te, pct_cyp_loci_with_te=round(100 * loci_te / len(loci), 1) if loci else 0,
               cyp_loci_with_any_hit=loci_any, te_per_cyp_locus=round(te_total / len(loci), 2) if loci else 0,
               te_per_10kb_window=round(10 * te_total / kb, 3) if kb else 0)
    for t in TOPS: row[f"te_hits_{t}"] = by_top.get(t, 0)
    row["te_hits_other"] = te_total - sum(by_top.get(t, 0) for t in TOPS)
    # xenobiotic targets: a target is present if it is one of a locus's Name values
    xp = xt = 0
    for t in targets:
        gl = [g for g, v in loci.items() if t.upper() in (n.upper() for n in v["names"])]
        te = sum(per_locus[g]["TE"] for g in gl); nt = sum(per_locus[g]["nonTE"] for g in gl)
        xp += bool(gl); xt += te
        xg_rows.append(dict(species=lab, exposure_group=row["exposure_group"], gene=t, present=bool(gl),
                            loci=len(gl), te_hits=te, simple_lowcomplex_hits=nt, has_te=te > 0))
    row.update(xeno_targets_present=xp, xeno_targets_with_te=sum(1 for r in xg_rows if r["species"] == lab and r["has_te"]), xeno_te_hits=xt)
    for g, v in loci.items():
        c = per_locus[g]
        cg_rows.append(dict(species=lab, exposure_group=row["exposure_group"], locus_id=g, names=",".join(v["names"]),
                            seq=v["seq"], start=v["start"], end=v["end"], te_hits=c["TE"], simple_lowcomplex_hits=c["nonTE"],
                            **{f"te_{t}": c.get("TE_" + t, 0) for t in TOPS}))
    sp_rows.append(row)

def write(name, rows):
    with open(f"{OUT}/{name}", "w", newline="") as o:
        w = csv.DictWriter(o, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
write("per-species.tsv", sp_rows); write("xenobiotic-genes.tsv", xg_rows); write("cyp-genes.tsv", cg_rows)
print(len(sp_rows), "species;", len(xg_rows), "xeno rows;", len(cg_rows), "cyp loci")
