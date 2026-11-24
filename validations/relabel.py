# NEW (2026-11-24): applies column_labels (LABELS, TE_KEEP) to an already-built cyp-te-validation.xlsx: the same
# result as a fresh make-validation-sheet.py run, without re-reading the RepeatMasker files.
import sys
from openpyxl import load_workbook
from openpyxl.styles import Font, Alignment
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from column_labels import LABELS, TE_KEEP
p = "/bioinformatics/atallah-repeats/summary/cyp-te-validation.xlsx"
wb = load_workbook(p)
inv = {v: k for k, v in LABELS.items()}
for t, keep in (("Cyp", None), ("Transposons", TE_KEEP)):
    ws, cw = wb[t], wb[f"{t} - columns"]
    short = [inv.get(c.value, c.value) for c in ws[1]]
    data = [dict(zip(short, r)) for r in ws.iter_rows(min_row=2, values_only=True)]
    doc = {}
    for r in cw.iter_rows(min_row=2, values_only=True):
        s = r[2] if len(r) > 2 and r[2] else inv.get(r[0], r[0]); doc[s] = r[1]
    cols = keep or short
    i = wb.sheetnames.index(t); del wb[t]; del wb[f"{t} - columns"]
    ws = wb.create_sheet(t, i); ws.append([LABELS.get(c, c) for c in cols])
    def num(v):
        if isinstance(v, str):
            try: return int(v)
            except ValueError:
                try: return float(v)
                except ValueError: return v
        return v
    for d in data: ws.append([num(d[c]) for c in cols])
    for c in ws[1]: c.font = Font(bold=True); c.alignment = Alignment(wrap_text=True, vertical="top")
    ws.freeze_panes = "B2"; ws.auto_filter.ref = ws.dimensions
    for k in range(1, len(cols) + 1): ws.column_dimensions[ws.cell(1, k).column_letter].width = 22
    ws.column_dimensions["A"].width = 34
    cw = wb.create_sheet(f"{t} - columns", i + 1); cw.append(["column", "what it is", "short name"])
    for c in cols: cw.append([LABELS.get(c, c), doc[c], c])
    for c in cw[1]: c.font = Font(bold=True)
    cw.column_dimensions["A"].width = 45; cw.column_dimensions["B"].width = 130; cw.column_dimensions["C"].width = 28
    for row in cw.iter_rows(min_row=2):
        for c in row: c.alignment = Alignment(wrap_text=True, vertical="top")
wb.save(p); print("relabelled", p, wb.sheetnames)
