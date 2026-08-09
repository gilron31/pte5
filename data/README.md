# Factor-base artifacts

One-time work (see `drafts/implementation_design.md` §3.1). Built once per
bound, then `mmap`'d and streamed; never resident.

## Format

Flat array of records, no header:

```
struct Record { uint32_t a; uint32_t b; };   // little-endian, 8 bytes
```

- One record per split prime `p ≡ 1 (mod 4)`, **ascending in `p`**.
- `a >= b > 0` and `p = a*a + b*b`. **The prime is not stored** — it is
  recoverable, so storing it would double the file for nothing.
- Entry count = filesize / 8.

Reading the whole base is therefore a `mmap` plus a cast; there is no decode
step, and the access pattern is pure sequential.

## Contents

| File | `P_max` | for bound `B = 65*P_max` | records | size |
| --- | --- | --- | --- | --- |
| `factorbase_1e12.bin` | 1.54e10 | 1e12 | 343,618,258 | 2.6 GB |

## Provenance

Segmented sieve + Tonelli–Shanks + Cornacchia, single-threaded, 2026-08-09:

```
primes to 1.54e10   687,242,934
split primes        343,618,258
segmented sieve            74.8 s
Cornacchia                462.1 s
total                     536.9 s   (8.9 min)
peak resident               0.2 MB
checksum              23649257516598   (sum of re^im over all records)
```

Verified after writing: size exact, first records (2,1),(3,2),(4,1),(5,2),(6,1)
→ p = 5,13,17,29,37; last p = 15,399,999,997; 400 random records all satisfy
`a>=b>0`, `p ≡ 1 (mod 4)`, `p` prime; ascending over a sampled run.

## Rebuilding

The builder is not yet in the tree — this was produced by a prototype pending
the segmented-sieve fix (**C8**) in `split_primes_up_to`. Until that lands,
treat this file as the input rather than something to regenerate on demand.
