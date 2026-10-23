"""N19: re-sieve the 3 quotients that the fixed MW sieve did not certify with the
default pool (primes < 120), now with primes < 220 and up to 24 steps.
Usage: ~/sage_venv/bin/python -u audit_retry_wide.py > retry3_wide_pool.log
"""
import sys, time, multiprocessing as mp
sys.path.insert(0, '/home/gilro/code/pte5/drafts/ag_research')
def one(args):
    label, rho, S, symes = args
    from sage.all import QQ
    from lift_pipeline import odd_model, chabauty, mw_sieve
    t0 = time.time()
    signed = [QQ(r) for r in rho] + [-QQ(r) for r in rho]
    e = [-signed[i - 1 if symes else i] for i in S]
    A, _ = odd_model(e)
    ok, p0, rep = chabauty(A, [QQ(0)])
    cert, surv, log, _ = mw_sieve(A, [QQ(0)], p0, rep, pool_max=220, max_steps=24)
    return f'{label} S={S}: wide-pool fixed sieve {"CERT" if cert else "NOT"} ({len(surv)} surv; {log[-1]}; {time.time()-t0:.0f} s)'
if __name__ == '__main__':
    from certify_parallel import tasks
    from chabauty_g2 import R_SYMES, CERTIFIED
    r = dict(tasks())
    todo = [('cat12781559905 point', r['cat12781559905 point'], [0, 6, 7, 9, 12, 15], False),
            ('F1k4 twin', r['F1k4 twin'], [0, 3, 6, 10, 13, 14], False),
            ('N15 Symes #2', R_SYMES, CERTIFIED[2], True)]
    with mp.get_context('fork').Pool(3) as pool:
        for line in pool.imap_unordered(one, todo):
            print(line, flush=True)
