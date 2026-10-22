"""A5, parallel: bottom-lift certificates for all remaining seeds (catalog to 1e12
and the Bremner family members in family_seeds.json), both twin forms, with a
pool of worker processes. Seeds already certified in certify_catalog.log or
certify_families.log are skipped.
Usage: ~/sage_venv/bin/python -u certify_parallel.py [workers] [tries]
"""
import json
import multiprocessing as mp
import re
import sys
import time


def done_labels():
    done = set()
    for path in ('certify_catalog.log', 'certify_families.log', 'certify_parallel.log'):
        try:
            for line in open(path):
                if 'CERTIFIED' not in line:
                    continue
                m = re.match(r'L1 =\s+(\d+)\s+(\w+):', line)
                if m:
                    done.add(f'cat{m.group(1)} {m.group(2)}')
                m = re.match(r'\s*(F\dk-?\d+)\s+(\w+)', line)
                if m:
                    done.add(f'{m.group(1)} {m.group(2)}')
                m = re.match(r'(\S+ \w+): CERTIFIED', line)
                if m:
                    done.add(m.group(1))
        except FileNotFoundError:
            pass
    return done


def tasks():
    from certify_catalog import seeds
    out = [(f'cat{L1} {form}', rho) for L1, form, rho in seeds('catalog_1e12.csv')]
    out += [(label, rho) for label, rho in json.load(open('family_seeds.json')).items()]
    return out


def work(args):
    label, rho, tries = args
    from lift_pipeline import certify_seed
    t0 = time.time()
    try:
        res = certify_seed(rho, tries=tries, verbose=False)
    except Exception as exc:
        return f'{label}: ERROR {type(exc).__name__}: {exc}'
    dt = time.time() - t0
    if res:
        return f'{label}: CERTIFIED no E_5 above (quotient {res[0]}, p0 = {res[1]}; {dt:.0f} s; max root {max(rho)})'
    return f'{label}: undecided after {tries} subsets ({dt:.0f} s; max root {max(rho)})'


if __name__ == '__main__':
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    tries = int(sys.argv[2]) if len(sys.argv) > 2 else 80
    done = done_labels()
    todo = [(label, rho, tries) for label, rho in tasks() if label not in done]
    print(f'{len(done)} seeds already certified; {len(todo)} to do with {workers} workers', flush=True)
    with mp.get_context('fork').Pool(workers) as pool:
        for line in pool.imap_unordered(work, todo):
            print(line, flush=True)
    print('all done', flush=True)
