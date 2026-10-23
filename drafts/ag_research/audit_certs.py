"""Audit of the G5 certificates (N17, N15), after the sieve fix of 2026-10-23.

For every logged certificate (seed form, quotient S):
  1. the seed is a nondegenerate E_4: prod (y - r_k^2) = (((y - L1)^2 - L2)^2 - L3)^2 - L4,
     distinct nonzero |r_k|;
  2. rank 1 re-proved (g2_rank1_generic.gp), Chabauty re-run (same p0?);
  3. the OLD sieve (no projection): do its primes have odd l | #J(F_q) with l not
     dividing ord(D0)? (the gap: the index of <D0> + J[2] was never tested there);
  4. the FIXED sieve (projection) on the same quotient: still certified?
  5. independent point search: hyperellratpoints on y^2 = prod (t - e_i), height H;
     every point found must be a known one (t = e_i);
  6. a second, independent certificate from a different random stream (seed = 2).
Usage: ~/sage_venv/bin/python -u audit_certs.py [workers] > audit_certs.log
"""
import json
import multiprocessing as mp
import re
import sys
import time

H = 3000


def logged():
    out = []
    for path in ('certify_catalog.log', 'certify_families.log', 'certify_parallel.log'):
        for line in open(path):
            if 'CERTIFIED' not in line:
                continue
            S = [int(v) for v in re.search(r'quotient \[([\d, ]+)\]', line).group(1).split(',')]
            p0 = int(re.search(r'p0 = (\d+)', line).group(1))
            m = re.match(r'L1 =\s+(\d+)\s+(\w+):', line)
            if m:
                label = f'cat{m.group(1)} {m.group(2)}'
            else:
                label = re.match(r'\s*(\S+ \w+)', line).group(1)
            out.append((label, S, p0))
    return out


def is_e4(rho):
    from sage.all import PolynomialRing, QQ
    R = PolynomialRing(QQ, 'y')
    y = R.gen()

    def psqrt(F):
        """monic G with deg(F - G^2) < deg G, and the remainder F - G^2"""
        d = F.degree() // 2
        G = y**d
        for k in range(d - 1, -1, -1):
            G += (F - G**2)[d + k] / 2 * y**k
        return G, F - G**2

    F = R.prod(y - QQ(r)**2 for r in rho)
    for _ in range(3):                  # E_4 -> E_3 -> E_2 -> E_1: each step F = G^2 - const
        F, rem = psqrt(F)
        if rem.degree() > 0:
            return False
    absr = [abs(r) for r in rho]
    return F.degree() == 1 and 0 not in absr and len(set(absr)) == len(absr)


def work(args):
    label, rho, S, p0_logged, symes = args
    from sage.all import QQ, pari
    from lift_pipeline import rank1, odd_model, chabauty, mw_sieve, certify_seed
    t0 = time.time()
    out = [f'{label} S={S}']
    try:
        out.append(f'E4 {"ok" if is_e4(rho) else "FAIL"}')
        signed = [QQ(r) for r in rho] + [-QQ(r) for r in rho]
        vals = [signed[i - 1 if symes else i] for i in S]
        e = [-v for v in vals]
        out.append(f'rank1 {"ok" if rank1(e) else "FAIL"}')
        A, to_X = odd_model(e)
        ok, p0, report = chabauty(A, [QQ(0)])
        out.append(f'chabauty {"ok" if ok else "FAIL"} p0={p0}' + ('' if symes or p0 == p0_logged else f' (logged {p0_logged})'))
        c_old, surv_old, log_old, (L_old, used_old) = mw_sieve(A, [QQ(0)], p0, report, project=False)
        from lift_pipeline import Jp, poly_from
        from sage.all import ZZ
        f = poly_from(A)
        W0 = QQ(f(0)).sqrt()
        gap = []
        for d in used_old:
            Jj = Jp(f, d['p'])
            o = Jj.order(Jj.cls(0, W0))
            k = 1
            for l, v in ZZ(Jj.N).factor():
                if l != 2 and o % l:
                    k *= l**v
            if k > 1:
                gap.append((d['p'], k))
        out.append(f'old sieve {"cert" if c_old else "NOT"} primes {[d["p"] for d in used_old]} untested-l parts {gap}')
        c_new, surv_new, log_new, (L_new, used_new) = mw_sieve(A, [QQ(0)], p0, report)
        out.append(f'fixed sieve {"CERT" if c_new else "NOT"} ({len(surv_new)} surv; primes {[d["p"] for d in used_new]}; '
                   + log_new[-1] + ')')
        pol = '*'.join(f'(x-({v}))' for v in e)
        pts = pari(f'hyperellratpoints({pol}, {H})')
        extra = [str(P[0]) for P in pts if QQ(P[0]) not in set(e)]
        out.append(f'ratpoints(H={H}) {len(pts)} found, extra {extra}')
        if not symes:
            res = certify_seed(rho, tries=80, seed=2, verbose=False)
            out.append(f'second cert {"CERT " + str(res[0]) + " p0=" + str(res[1]) if res else "NONE"}')
    except Exception as exc:
        out.append(f'ERROR {type(exc).__name__}: {exc}')
    out.append(f'{time.time() - t0:.0f} s')
    return ' | '.join(out)


if __name__ == '__main__':
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    from certify_parallel import tasks
    rhos = dict(tasks())
    todo = [(label, rhos[label], S, p0, False) for label, S, p0 in logged()]
    from chabauty_g2 import R_SYMES, CERTIFIED
    todo += [(f'N15 Symes #{i}', R_SYMES, S, None, True) for i, S in enumerate(CERTIFIED)]
    print(f'{len(todo)} certificates to audit', flush=True)
    with mp.get_context('fork').Pool(workers) as pool:
        for line in pool.imap_unordered(work, todo):
            print(line, flush=True)
    print('audit done', flush=True)
