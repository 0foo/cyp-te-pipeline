# NEW (not a lab file, 2026-11-22): pairs each stage 1 RepeatMasker sample with its renamed annotation,
# and checks that the two describe the same assembly. It writes pairs.tsv, which run-pairing.sh,
# make-exposure.py and make-summary.py read.
#
# WHY: Locate_TE.py compares seqids and coordinates between the annotation and the .rm.out. If the
#   annotation belongs to a different assembly of the species, nothing matches, or worse, the wrong
#   coordinates match, and the run still finishes without any error.
# HOW
#   label  the sample name normalised (accession or 'nanopore' suffix dropped, non-alphanumerics -> _,
#          upper case), matched to the longest renamed-GFF label it equals or starts with.
#   checks the sequence lengths implied by the .rm.out (end + left); the fraction of the GFF's genes on
#          sequences present in the .rm.out; genes ending past a sequence's end; and Cyp genes (any
#          Reg_Gene_Full.txt name in the attributes, the same substring rule as repeatOpp.py) and how
#          many of them are on .rm.out sequences.
#   status RUN (fraction >= 0.90 and nothing past an end) | SKIP_WRONG_ASSEMBLY (fraction < 0.50: the
#          annotation is of a different assembly) | CHECK (in between) | NO_ANNOTATION |
#          DUPLICATE_LABEL (more than one sample would RUN for one annotation; choose by hand).
#   The full run gave 301 RUN, 3 SKIP_WRONG_ASSEMBLY and 15 NO_ANNOTATION.
# INPUTS: 3-pair-te-with-cyp/input/rm-out/*.rm.out and input/renamed-gff/*_final_withDmelNames.gff;
#   D. melanogaster uses its own un-renamed GFF.
# OUTPUT: 3-pair-te-with-cyp/pairs.tsv (rm_sample, gff_label, gene_frac_on_rm_contigs,
#   genes_past_seq_end, cyp_genes, cyp_on_rm_contigs, status), plus a summary on stdout.
# Usage: python3 -B make-pairs.py
import os, re, csv, collections
A = "/bioinformatics/atallah-repeats"
H = f"{A}/3-pair-te-with-cyp"
RM, GD = f"{H}/input/rm-out", f"{H}/input/renamed-gff"
CYP = [l.strip() for l in open("/root/projects/cyp-te-pipeline/data-analysis/1-pair-te-with-cyp/Reg_Gene_Full.txt")
       if l.strip() and not l.strip().startswith("D")]
gff_path = {f[:-len("_final_withDmelNames.gff")]: f"{GD}/{f}" for f in os.listdir(GD)}
gff_path["DROSOPHILA_MELANOGASTER"] = f"{A}/2-rename-genes/input/gffs/DROSOPHILA_MELANOGASTER_final.gff"

def norm(s):
    s = re.split(r"\.(GC[AF]_|nanopore)", s)[0]
    return re.sub(r"_+", "_", re.sub(r"[^A-Za-z0-9]+", "_", s)).strip("_").upper()

def label_for(s):
    n = norm(s)
    c = [g for g in gff_path if n == g or n.startswith(g + "_")]
    return max(c, key=len) if c else None

def seq_lengths(p):
    L = {}
    with open(p) as f:
        for i, line in enumerate(f):
            x = line.split()
            if i >= 3 and len(x) >= 8:
                L[x[4]] = int(x[6]) + int(x[7].strip("()"))
    return L

rows = []
for f in sorted(os.listdir(RM)):
    s = f[:-len(".rm.out")]
    g = label_for(s)
    if not g:
        rows.append(dict(rm_sample=s, gff_label="-", status="NO_ANNOTATION")); continue
    L = seq_lengths(f"{RM}/{f}")
    n = on = past = nc = cyp_on = 0
    with open(gff_path[g]) as fh:
        for line in fh:
            x = line.split("\t")
            if len(x) > 8 and x[2] == "gene":
                n += 1; ok = x[0] in L; on += ok
                past += ok and int(x[4]) > L[x[0]]
                if any(c in x[8] for c in CYP): nc += 1; cyp_on += ok
    frac = on / n if n else 0
    st = "RUN" if frac >= 0.90 and not past else "SKIP_WRONG_ASSEMBLY" if frac < 0.50 else "CHECK"
    rows.append(dict(rm_sample=s, gff_label=g, gene_frac_on_rm_contigs=round(frac, 4), genes_past_seq_end=past,
                     cyp_genes=nc, cyp_on_rm_contigs=cyp_on, status=st))
runs = collections.Counter(r["gff_label"] for r in rows if r["status"] == "RUN")
for r in rows:
    if r["status"] == "RUN" and runs[r["gff_label"]] > 1: r["status"] = "DUPLICATE_LABEL"
cols = ["rm_sample", "gff_label", "gene_frac_on_rm_contigs", "genes_past_seq_end", "cyp_genes", "cyp_on_rm_contigs", "status"]
with open(f"{H}/pairs.tsv", "w", newline="") as o:
    w = csv.DictWriter(o, fieldnames=cols, delimiter="\t", restval="", lineterminator="\n"); w.writeheader(); w.writerows(rows)
print(collections.Counter(r["status"] for r in rows))
for r in rows:
    if r["status"] not in ("RUN", "NO_ANNOTATION"): print({k: r.get(k) for k in cols})
used = {r["gff_label"] for r in rows if r["status"] == "RUN"}
print("annotations with no RUN sample:", sorted(set(gff_path) - used))
