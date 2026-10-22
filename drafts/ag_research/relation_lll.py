"""Q1: short integer linear relations among the 8 signed roots, via LLL.

For a random 8-vector of size N the relation lattice has rank 7 and covolume
about N, so its successive minima are about N^(1/7). An anomalously short
relation (besides the slice relation A - C - E + G = 0) would suggest a
family-1-type hyperplane section through the point. Family members are the
controls. Roots are taken in the signed slice labeling (A..H) from invert.py.
"""
import subprocess

from known_e4 import all_known, twin
from invert import signed_labelings


def short_relations(v):
    gp = f'''
v = {list(v)};
K = matkerint(Mat(v));
L = K * qflll(K);
print(vector(#L, j, sqrtint(norml2(L[, j]))));
for (j = 1, 3, print(L[, j]~));
'''
    return subprocess.run(['gp', '-q'], input=gp, capture_output=True, text=True).stdout.strip()


if __name__ == '__main__':
    rows = all_known(kmax=7)
    picks = ['cat33065996717', 'F1k1', 'F1k2', 'F1k3', 'F2k6', 'F2k7', 'F2k8', 'cat12781559905', 'cat17247059713']
    for name in picks:
        r = rows[name]
        for form, rr in (('root', r), ('twin', twin(r))):
            labs = list(signed_labelings(rr))
            if labs:
                v = labs[0]
                n = max(abs(c) for c in v)
                print(f'{name} ({form} form, on slice), max root {n}, N^(1/7) = {n ** (1/7):.1f}')
                print('   ' + short_relations(v).replace('\n', '\n   '))
                break
        else:
            v = r
            n = max(abs(c) for c in v)
            print(f'{name} (off slice), max root {n}, N^(1/7) = {n ** (1/7):.1f}')
            print('   ' + short_relations(v).replace('\n', '\n   '))
