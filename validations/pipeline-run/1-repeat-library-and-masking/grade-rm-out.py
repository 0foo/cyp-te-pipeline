# NEW (not a lab file, 2026-11-22): grades the JSON written by check-rm-out.py and gives each
# RepeatMasker sample a verdict of OK or BAD.
#
# A sample is BAD if any of the following holds:
#   - the .out file is empty, contains binary garbage or NUL bytes, or has unparseable rows (up to 20
#     "repeatpos" rows are tolerated, because RepeatMasker itself writes a few negative positions
#     such as "(-2)")
#   - the header is bad, or the file has no final newline (truncated)
#   - rows disagree on a sequence's length, or there are no data rows
#   - the union of hit intervals differs from the .tbl "bases masked" by more than 2%
#   - the .rm.tbl is missing, the masker log does not end in "done", or it has error-like lines
# OUTPUT: TSV on stdout (sample, verdict, reasons) and a {verdict: count} summary on stderr.
# Usage: python3 -B grade-rm-out.py RESULT.json > grade.tsv
import sys, json, collections

res = json.load(open(sys.argv[1]))
counts = collections.Counter()
print('sample\tverdict\treasons')
for r in res:
    why = []
    if r['status'] == 'EMPTY':
        why.append('empty file')
    else:
        if r.get('badblocks'): why.append(f"binary garbage in {r['badblocks']} blocks from byte {r['first_bad']}")
        if r.get('nul'): why.append(f"{r['nul']} NUL bytes")
        bad = {k: v for k, v in r.get('bad', {}).items() if k != 'garbage'}
        # RepeatMasker itself emits a few rows with negative repeat positions, e.g. "(-2)"; tolerate up to 20
        if bad.get('repeatpos', 0) <= 20: bad.pop('repeatpos', None)
        if bad: why.append('unparseable rows ' + ','.join(f'{k}={v}' for k, v in bad.items()))
        if not r.get('header_ok'): why.append('bad header')
        if not r.get('ends_nl'): why.append('truncated (no final newline)')
        if r.get('implied_mismatch'): why.append(f"{r['implied_mismatch']} rows disagree on sequence length")
        if r.get('data', 0) == 0: why.append('no data rows')
        # masked bases in .tbl vs union of hit intervals in .out (they should agree closely)
        m, c = r.get('tbl_masked'), r.get('cover')
        if m and c and abs(m - c) / m > 0.02: why.append(f'coverage {c} vs tbl masked {m}')
    if not r.get('tbl_masked'): why.append('rm.tbl missing/unreadable')
    if not r.get('mlog_done'): why.append(f"masker log does not end in done ({r.get('mlog_last', r.get('mlog'))!r})")
    if r.get('mlog_nbad'): why.append(f"masker log has {r['mlog_nbad']} error-like lines")
    v = 'OK' if not why else 'BAD'
    counts[v] += 1
    print(f"{r['sample']}\t{v}\t{'; '.join(why)}")
print(dict(counts), file=sys.stderr)
