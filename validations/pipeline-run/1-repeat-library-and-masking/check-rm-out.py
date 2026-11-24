# NEW (not a lab file, 2026-11-22): integrity check of RepeatMasker output.
# Usage: python3 check-rm-out.py RESULT.json [RMODELER_ROOT]   (root holds out/ logs/ state/)
import sys, os, re, json
from multiprocessing import Pool

ROOT = sys.argv[2] if len(sys.argv) > 2 else '/bioinformatics/atallah-repeats/rmodeler'
OUT = ROOT + '/out'
LOGS = ROOT + '/logs'
STATE = ROOT + '/state'
LEFT = re.compile(r'^\((\d+)\)$')
INT = re.compile(r'^\d+$')
FLT = re.compile(r'^\d+(\.\d+)?$')
NONTE = ('Simple_repeat', 'Low_complexity', 'Satellite', 'rRNA', 'tRNA', 'snRNA', 'scRNA', 'srpRNA', 'ARTEFACT')


def parse_tbl(p):
    r = {'tbl_size': os.path.getsize(p)}
    try:
        t = open(p, errors='replace').read()
    except Exception as e:
        r['tbl_err'] = str(e); return r
    m = re.search(r'sequences:\s+(\d+)', t); r['tbl_nseq'] = int(m.group(1)) if m else None
    m = re.search(r'total length:\s+(\d+) bp', t); r['tbl_len'] = int(m.group(1)) if m else None
    m = re.search(r'bases masked:\s+(\d+) bp', t); r['tbl_masked'] = int(m.group(1)) if m else None
    m = re.search(r'file name:\s+(\S+)', t); r['tbl_file'] = m.group(1) if m else None
    return r


def parse_log(s):
    r = {}
    p = f'{LOGS}/{s}.masker.log'
    if not os.path.exists(p):
        r['mlog'] = 'missing'; return r
    t = open(p, errors='replace').read()
    r['mlog_size'] = len(t)
    tail = t.rstrip().splitlines()[-1] if t.strip() else ''
    r['mlog_last'] = tail[:80]
    r['mlog_done'] = tail.strip() == 'done'
    bad = [l for l in t.splitlines() if re.search(r'error|kill|die[sd]?\b|terminat|abort|cannot|can\'t|failed|warning|out of memory|segmentation', l, re.I)]
    r['mlog_bad'] = bad[:5]; r['mlog_nbad'] = len(bad)
    m = re.findall(r'in batch (\d+) of (\d+)', t)
    if m:
        tot = int(m[0][1]); seen = {int(a) for a, b in m}
        r['mlog_batches'] = tot; r['mlog_batches_seen'] = len(seen)
    r['mlog_has_tbl_stage'] = 'Generating table output' in t
    p2 = f'{LOGS}/{s}.log'
    if os.path.exists(p2):
        t2 = open(p2, errors='replace').read()
        bad2 = [l for l in t2.splitlines() if re.search(r'error|killed|terminat|abort|out of memory|segmentation|traceback', l, re.I)]
        r['log_bad'] = bad2[:5]; r['log_nbad'] = len(bad2)
    else:
        r['log'] = 'missing'
    for k in ('done', 'masked', 'failed'):
        r['state_' + k] = os.path.exists(f'{STATE}/{k}/{s}')
    return r


BADB = re.compile(rb'[^\x09\x0a\x20-\x7e]')
import math, collections
def blockscan(p, B=4096):
    runs = []; nbad = 0; i = 0; first = None; ent = None
    with open(p, 'rb') as f:
        while True:
            b = f.read(B)
            if not b: break
            isbad = bool(BADB.search(b))
            if isbad:
                nbad += 1
                if first is None:
                    m = BADB.search(b); first = i * B + m.start()
                    c = collections.Counter(b); ent = -sum(v/len(b)*math.log2(v/len(b)) for v in c.values())
                if runs and runs[-1][1] == i: runs[-1][1] = i + 1
                else: runs.append([i, i + 1])
            i += 1
    return dict(nblocks=i, badblocks=nbad, first_bad=first, bad_runs=len(runs),
                runs_mib=[(round(a*B/1048576, 3), round(b*B/1048576, 3)) for a, b in runs[:6]], ent_first=ent)

def check(s):
    p = f'{OUT}/{s}.rm.out'
    r = {'sample': s, 'size': os.path.getsize(p)}
    r.update(parse_tbl(f'{OUT}/{s}.rm.tbl'))
    for ext in ('-families.fa', '-families.stk'):
        q = f'{OUT}/{s}{ext}'
        r['sz' + ext] = os.path.getsize(q) if os.path.exists(q) else None
    r.update(parse_log(s))
    if r['size'] == 0:
        r['status'] = 'EMPTY'; return r
    r.update(blockscan(p))
    nul = nonascii = 0
    bad = {}; badex = []
    nlines = ndata = 0
    seqlen = {}; seqorder = []; last_seq = None; reappear = 0
    implied_mismatch = 0; impl_ex = []
    cover = 0; cur_s = None; cur_a = cur_b = None
    unsorted = 0
    classes = {}
    classified = 0
    lastid = 0; id_back = 0
    with open(p, 'rb') as f:
        data_iter = f
        header = []
        ends_nl = True
        for raw in data_iter:
            nlines += 1
            if b'\x00' in raw: nul += raw.count(b'\x00')
            try:
                line = raw.decode('ascii')
            except UnicodeDecodeError:
                nonascii += 1; line = raw.decode('latin-1')
            ends_nl = raw.endswith(b'\n')
            if BADB.search(raw):
                bad['garbage'] = bad.get('garbage', 0) + 1
                continue
            if nlines <= 3:
                header.append(line.rstrip('\n')); continue
            if not line.strip():
                continue
            fl = line.split()
            why = None
            if not (14 <= len(fl) <= 16): why = f'nfields={len(fl)}'
            else:
                if not INT.match(fl[0]): why = 'score'
                elif not all(FLT.match(x) for x in fl[1:4]): why = 'div/del/ins'
                elif not (INT.match(fl[5]) and INT.match(fl[6])): why = 'begin/end'
                elif not LEFT.match(fl[7]): why = 'left'
                elif fl[8] not in ('+', 'C'): why = 'strand'
                elif not fl[10] or fl[10].isdigit(): why = 'class'
                else:
                    b, e = int(fl[5]), int(fl[6])
                    if b > e or b < 1: why = 'begin>end'
                    else:
                        rp = fl[11:14]
                        if fl[8] == '+':
                            ok = INT.match(rp[0]) and INT.match(rp[1]) and LEFT.match(rp[2])
                        else:
                            ok = LEFT.match(rp[0]) and INT.match(rp[1]) and INT.match(rp[2])
                        if not ok: why = 'repeatpos'
                        elif len(fl) >= 15 and not INT.match(fl[14]): why = 'ID'
                        elif len(fl) == 16 and fl[15] != '*': why = 'col16'
            if why:
                bad[why] = bad.get(why, 0) + 1
                if len(badex) < 3: badex.append(f'L{nlines}: {line.rstrip()[:160]}')
                continue
            ndata += 1
            q = fl[4]; b = int(fl[5]); e = int(fl[6]); left = int(fl[7][1:-1])
            L = e + left
            if q != last_seq:
                if q in seqlen: reappear += 1
                else: seqorder.append(q)
                # flush coverage
                if cur_a is not None: cover += cur_b - cur_a + 1
                cur_s = q; cur_a = cur_b = None
                last_seq = q
            if q in seqlen:
                if seqlen[q] != L:
                    implied_mismatch += 1
                    if len(impl_ex) < 3: impl_ex.append(f'L{nlines} {q} {seqlen[q]} vs {L}')
            else:
                seqlen[q] = L
            if cur_a is None:
                cur_a, cur_b = b, e
            else:
                if b < cur_a: unsorted += 1
                if b <= cur_b + 1:
                    cur_b = max(cur_b, e)
                else:
                    cover += cur_b - cur_a + 1; cur_a, cur_b = b, e
            c = fl[10]
            top = c.split('/')[0]
            classes[top] = classes.get(top, 0) + 1
            if not (top in NONTE or top.startswith('Unknown') or top == 'Unspecified'):
                classified += 1
            if len(fl) >= 15:
                i = int(fl[14])
                if i < lastid - 0: id_back += 1
                lastid = max(lastid, i)
    if cur_a is not None: cover += cur_b - cur_a + 1
    r.update(dict(lines=nlines, data=ndata, header=header, header_ok=(len(header) == 3 and 'SW' in header[0] and 'score' in header[1] and header[2].strip() == ''),
                  ends_nl=ends_nl, nul=nul, nonascii_lines=nonascii, bad=bad, badex=badex,
                  nseq_out=len(seqlen), sum_implied=sum(seqlen.values()), max_implied=max(seqlen.values()) if seqlen else 0,
                  implied_mismatch=implied_mismatch, impl_ex=impl_ex, reappear=reappear, unsorted=unsorted,
                  cover=cover, classified=classified, classes=dict(sorted(classes.items(), key=lambda x: -x[1])[:8]),
                  maxid=lastid, id_back=id_back, last_line_seq=last_seq, last_seq_index=len(seqorder)))
    r['status'] = 'parsed'
    return r


if __name__ == '__main__':
    samples = sorted(f[:-7] for f in os.listdir(OUT) if f.endswith('.rm.out'))
    with Pool(8) as pool:
        res = pool.map(check, samples, chunksize=1)
    json.dump(res, open(sys.argv[1], 'w'), indent=1)
    print(len(res))
