"""Run chabauty_g2.chabauty on one certified quotient and print every disc."""
import sys
from chabauty_g2 import chabauty, CERTIFIED
i = int(sys.argv[1]) if len(sys.argv) > 1 else 0
p = int(sys.argv[2]) if len(sys.argv) > 2 else None
ok, p, report = chabauty(CERTIFIED[i], p=p, verbose=False)
print('S =', CERTIFIED[i], 'p =', p)
for Qbar, known, bound in report:
    print(f'   disc {Qbar}: known {known}, Strassmann bound {bound}')
print('total known', sum(r[1] for r in report), ' total bound', sum(r[2] for r in report))
