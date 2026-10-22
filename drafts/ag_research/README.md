# `drafts/ag_research/` — scripts behind `../ag_research_notes.md`

Each script implements one notebook entry (N1–N18). Run everything from this
directory. Toolchain:

- Python 3 with sympy and numpy;
- Singular 4.4;
- PARI/GP 2.17;
- Sage (passagemath 10.8.13) via `~/sage_venv/bin/python`;
- g++ with OpenMP, and nvcc 12.9 (see the build notes below).

**Shared modules.**

- `slice.py`: Bremner's slice in $(x:y:z)$, with the family curves and the
  degenerate factors.
- `known_e4.py`: all known $E_4$ (the catalog, Bremner's Table 1, the families),
  twin forms, identity keys.
- `phi_degree.py`: the invariants $L_1..L_4$ on the slice in $(s:t:u)$, and the
  $\varphi$-fiber script.
- `lift_pipeline.py`: the generic G5 certificate pipeline.

## Index

| entry | what | scripts | time |
| --- | --- | --- | --- |
| N1 | catalog to $10^{12}$ (`catalog_1e12.csv`, 15 $E_4$) | `code/native` `enum_campaign` | 558 s |
| N2–N3 | identities; small linear relations | `census.py` | 10 s |
| N4 | the exemplar on the slice; $(x:y:z)$ preimages | `invert.py` | 6 s |
| N4a–b | Q1 ansatz (family-2-type cubics), all 6 split lines | `ansatz_search.py`, `ansatz_exemplar_all.sh` | 15 min |
| N5, N10 | slice points to height 170 / 500 / 1200 (Python) | `slice_search.py`, `slice_search_big.py`, `slice_sieve.py`, `degen_probe.py`, `degen_curves.py` | 3 s / 21 s / ~5 min |
| N6 | uniqueness of Bremner's plane | `special_planes.py` | 1 s |
| N7 | top-node local ceilings, $p \le 53$ | `top_node_witness.py` (uses `../ag_pilots/local_upper.py`) | 1 s |
| N8 | family curves in $(s:t:u)$; family 1's image in $(J_1, J_2)$ | `fam_stu.sing`, `implicit_fam1_stu.sing` | 1 s |
| N8a | no $E_5$ inside family 1 | `sing_image.sing`, `sing_rational_fam1.gp`, `twin_check.py fam1 ...` | 1 s |
| N8c | no $E_5$ across the families | `cross_family.gp` (7 min), `cross_images.gp` (0.1 s) | |
| N8d | family 2's image (modular, exact check); no $E_5$ inside family 2 | `modular_image.py fam2`, `verify_image.py`, `exact_check_fam2.gp`, `sing_rational_fam2.gp`, `twin_check.py fam2 ...` | 18 s + 33 s |
| N8b | lines / conics tangent to the branch conics | `tangent_lines.py`, `bitangent_conics.sing` | 1 s |
| N9, N18 | $\varphi$ on the slice: degree, fibers, scans | `phi_degree.py`, `phi_fiber_pts.py`, `phi_scan_check.gp`, `phi_scan2.gp` (21 s), `phi_scan400.gp` (209 s), `member_partners.gp` | |
| N11 | LLL short relations | `relation_lll.py` | 1 s |
| N12 | $\mathrm{Sing}(V_5)$ strata; local types (slow) | `v5_singular.py` (0.6 s), `v5_local_type.py` (hours) | |
| N13 | genus-2 descent; validation; rank-1 certificates | `g2descent.gp`, `g2descent_test.gp`, `g2descent_test_all.gp`, `g2_lift_ranks.gp`, `g2_rank1_certify.gp`, `g2_rank1_generic.gp` | seconds |
| N14–N15 | Chabauty + MW sieve on Symes' quotients | `chabauty_smoke.py`, `chabauty_g2.py`, `chabauty_one.py`, `mw_sieve.py i` | 30 s per curve |
| N16 | calibration ranks; omission soundness test | `calib_ranks.gp`, `omission_test.py` | 25 s |
| N17 | certificates over the catalog / family members | `certify_catalog.py`, `certify_families.py` (+ `family_seeds.json`) | 20–80 s per seed |
| N18 | slice search to height 5000+ (C++ / CUDA) | `slice_search.cpp`, `slice_search_gpu.cu`, `classify_hits.py` | 1 s at $H = 500$; GPU 90 s at $H = 5000$ |

Files ending in `_tmp.sing`, `phi_fiber*.sing` and `modimg_tmp.sing` are scripts
generated on the fly and can be deleted. Logs (`*.log`, `hits_*.txt`) are run
outputs.

## Builds

```sh
g++ -O3 -march=native -fopenmp -o slice_search slice_search.cpp
nvcc -O3 -std=c++20 -ccbin=/usr/bin/g++-13 -arch=sm_89 -U_GNU_SOURCE -o slice_search_gpu slice_search_gpu.cu
```

`-U_GNU_SOURCE` is required on this machine. glibc 2.43 declares
`cospi`/`sinpi` `noexcept`, which clashes with CUDA 12.9's `math_functions.h`, so
even an empty `.cu` fails without it. `code/native`'s CUDA preset will need the
same flag when it is rebuilt.

## Operational rules (learned the hard way)

- **Do not use Singular `modStd`.** It forks a worker per core at gigabytes each,
  and `setcores` did not cap it. Use the sequential multi-modular route
  (`modular_image.py`) instead.
- **In PARI, set `parisizemax` only in the top-level script.** Inside a file
  passed to `read()` it aborts the rest of the read.
- **Never `pkill -f` with a pattern that appears in your own command line.** Kill
  by PID.
- **Singular reads `u^8/39304` as `u^(8/39304)`.** Clear denominators before
  generating Singular input.
- **Singular prints `3x2y` for $3x^2y$ under the default `short` setting.** Set
  `short = 0;` before parsing its output.
