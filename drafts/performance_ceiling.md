# Performance Ceilings — measured

**What this is.** `drafts/gpu_implementation_plan.md` (rev. 3) ranks the GPU's
value from cost *models*. Two of those models moved by three orders of magnitude
while the document was being written, which is a poor basis for building kernels.
This report replaces the models with measurements: where the CPU time actually
goes, how far each phase is from its hardware floor, and what the eight GPU
primitives the next milestones need can actually do on the device we have.

**The organising question is headroom, not speed.** "42 ns/norm" decides nothing;
"60% of the load-port ceiling" decides whether to keep working on it.

Everything here is reproducible from `code/native/bench/`. Run:

```sh
cmake --preset dev && cmake --build --preset dev
./build/dev/bench/bench_walk_phases 1000000000 3
./build/dev/bench/bench_meet 10000000 7
./build/dev/bench/gpu/bench_gpu_primitives
```

Machine: AMD Ryzen 7 8845HS (Zen 4, 8C/16T, 0.4–5.1 GHz) + RTX 4060 Laptop
(sm_89, 24 SMs, 8 GB, 128-bit). `kernel.perf_event_paranoid=1` and the
`performance` governor are set; `ncu` is still blocked (see §5).

---

## 0. Methodology, and why it needed building

Three facts made every prior number in this project soft:

1. **Wall clock varies up to 2×** on this laptop. One published claim
   (`implementation_design.md` §9.4) had to be retracted because of it.
2. **`perf` was blocked** at `perf_event_paranoid=4`.
3. **The GPU idles at 210 MHz of 3105.** An un-warmed kernel measurement is
   ~15× pessimistic — and wrong *reproducibly*, so repeating it confirms it.

`bench/harness.h` addresses all three:

- **Interleaving.** Variants run round-robin, never back-to-back, so a thermal
  ramp hits every variant equally and the *ratio* survives even when absolutes
  drift.
- **Counters, not clocks.** `bench/perf_counters.h` wraps `perf_event_open(2)`
  directly — scoped to a *region*, not a process, so per-phase IPC falls out of
  the same ablation that gives per-phase time. **cycles/unit is invariant under
  frequency scaling**, which removes the variance problem rather than working
  around it.
- **A dynamometer.** A dependent multiply-add chain of known cycle cost runs
  between rounds and reports the clock the sweep actually ran at. It reads
  5.0–5.1 GHz against a 5137 MHz hardware maximum — independently confirming
  it is correct — and stamps UNSTABLE on the table when drift exceeds 5%.
- **GPU clock discipline.** `bench/gpu/gpu_harness.cuh` warms until the SM clock
  plateaus and prints the achieved clock beside every result. This is not
  optional bookkeeping: P3 shows the clock *sagging to 1.09–1.37 GHz* under
  heavy sort traffic, which would otherwise have been silently attributed to the
  algorithm.

Two mistakes this discipline caught in its own benchmarks, both worth recording
because they are the standard ones:

- The dynamometer first read a **31% clock drift that never happened** — the
  0.4 ms window is short enough that one scheduler tick inflates wall time while
  the user-space cycle counter does not move. Fixed by measuring against thread
  CPU time and taking best-of-three (interference is one-sided).
- The GPU ALU benchmark first reported **1231% of the hardware roofline**,
  because `a += constant` in a loop is strength-reduced to `a += constant*reps`
  and the loop is deleted. Fixed with `asm volatile` PTX. *Any ALU measurement
  above 100% of roofline is a deleted loop, not a discovery.*

---

## 1. Where the CPU time goes (B1, ablation + counters)

Four instantiations of the same walk, each adding one phase; consecutive
differences are the phases. Self-checking: all four levels must report the same
norm count, and those counts are the 11,083 / 131,437 / 1,464,114 fixtures. The
skeleton row also reproduces §8's independently measured "DFS walk only" to
within 0.1 ns/norm.

| $B$ | DFS skeleton | `Spine::push` | ν extraction | `meet` | total ns/norm |
| --- | --- | --- | --- | --- | --- |
| $10^6$ | 12.6 (14.4%) | **44.6 (51.0%)** | 2.5 (2.9%) | 27.7 (31.7%) | 87.4 |
| $10^7$ | 11.8 (11.9%) | **41.0 (41.5%)** | 4.1 (4.2%) | 41.8 (42.4%) | 98.7 |
| $10^8$ | 10.9 (9.7%) | 41.1 (36.5%) | 4.1 (3.7%) | **56.5 (50.1%)** | 112.7 |
| $10^9$ | 10.4 (8.1%) | 41.2 (32.2%) | 3.9 (3.1%) | **72.5 (56.7%)** | 128.1 |

**The two halves behave completely differently.** `Spine::push` is *flat* —
44.6 → 41.2 ns/norm across three decades, an $O(1)$ per-norm cost. `meet` grows
27.7 → 72.5 and **crosses `push` between $10^7$ and $10^8$**.

That is the confirmation §2 of the design doc needed: the profile shape tracks
the $D>8$ work share (35.3% → 59.5% → 74.2%), so the $20.5\times$/decade growth
of meet work — and the $10^{12}$ extrapolation resting on it — is now confirmed
on real profiles rather than on quad-partition counts alone.

**Consequence: the optimization target has moved.** At the design bound the
campaign is a meet. Branch misses per norm in the meet phase grow 2.28 → 4.07 →
5.93 across $10^7$–$10^9$ while IPC falls 3.35 → 2.94; the meet is the only phase
whose behaviour degrades with scale.

---

## 2. CPU ceilings, per phase

### 2.1 The meet (B3) — load-bound at ~60% of ceiling

| | measured | ceiling | achieved |
| --- | --- | --- | --- |
| shipped `meet<>` dispatch | 180.6 cyc/norm | 109.9 | **60.9%** |
| all-hashed | 245.4 | — | 44.8% |
| all-direct ($D \le 16$) | 215.8 | — | 50.9% |

The ceiling model matters more than the number, and a first draft of it was
wrong in an instructive way. It is tempting to call `meet_direct` ALU-bound at
1 add + 1 compare per quad-partition with `sq[]` in registers. But
`partitions<D>()` is a constexpr **array**, not unrolled code, so each iteration
loads one 4-byte `Partition` and then indexes `sq[]` four times with runtime
indices — **five loads against three load ports (1.67 cyc), versus three ALU ops
against four pipes (0.75 cyc)**.

> **The direct meet is load-bound, by more than 2×.** The way past its ceiling is
> not wider adds; it is turning those indices from *data* into *code* — compile-
> time unrolling, or SIMD shuffles over a fixed permutation — so the four gathers
> disappear. That is a specific, testable proposal rather than "try SIMD".

Per-D, on the real corpus at $B = 10^7$ (`%ceil` against the load-bound model):

| $D$ | norms | direct ns | hashed ns | winner | %ceil direct |
| --- | --- | --- | --- | --- | --- |
| 4 | 101,605 | 6.0 | 19.0 | direct | 31.8% |
| 6 | 14,608 | 36.8 | 51.5 | direct | 46.6% |
| 8 | 11,305 | 149.3 | 148.8 | **tie** | 50.8% |
| 10 | 362 | 432.2 | 209.5 | hashed | 51.9% |
| 16 | 261 | 4012.5 | 1422.6 | hashed | 48.0% |

The shipped cut at $D \le 8$ is confirmed and is genuinely a coin flip at 8
(149.3 vs 148.8 ns), exactly as `implementation_design.md` §7.1 said the sweep
could not resolve.

### 2.2 `meet_hashed` — 82% of it is the sort

Decomposed on the $D>8$ tail (2,790 norms, 65 pair keys/norm). The staged
transcription is gated against the shipped function: hit counts must agree before
any number prints.

| stage | ns/norm | cyc/norm | share |
| --- | --- | --- | --- |
| hoist `sq[]` | 7.4 | 36.7 | 1.1% |
| key generation | 47.2 | 238.0 | 7.1% |
| **sort** | **542.8** | **2726.0** | **81.9%** |
| scan + emit | 65.4 | 326.5 | 9.9% |

And the sort carries **86 branch misses per norm** against 0.34 for key
generation, with IPC collapsing from 5.55 to 1.53. At ~13 cycles per mispredict
that is ~1120 of the sort's 2726 cycles — **41% of `meet_hashed`'s dominant stage
is branch misprediction.**

> **This is the highest-value open optimization in the project.** $D>8$ is 3% of
> norms and 74% of direct-meet work at $10^8$, `meet_hashed` is where it lands,
> and 82% of `meet_hashed` is a comparison sort that is mispredict-bound. Any
> sort-avoiding scheme — radix, or a small hash table finding equal keys directly
> — attacks the largest single cost at the design bound. Note this is the same
> structural idea as the value-windowing that made S-E 1000× cheaper
> (`gpu_implementation_plan.md` §3.1), which answers that document's open
> question 2 in the affirmative.

### 2.3 `Spine::push` — the sorting-network idea does *not* pay (B4)

`implementation_design.md` §10.1 found ~85% of `push` is `std::sort` on a
**three-element** array (mean 3.16 values, 90.3% are $n \le 4$), and §10.3
speculated ~10× headroom from a branch-free sorting network.

**That was wrong, and the counters say why.**

| variant | cyc/push | net of copy floor | branch misses/push |
| --- | --- | --- | --- |
| `std::sort` + `unique` (shipped) | 47.8 | 33.8 | **0.227** |
| branchless insertion + `unique` | 48.8 | 34.8 | 0.062 |
| sorting network + `unique` | 37.6 | **23.6** | 0.091 |
| copy only (bench overhead) | 14.0 | — | 0.013 |

Networks cover 99.1% of pushes and were verified to agree with `std::sort` on
every push in the corpus. They win **30% of the sort itself, ~21% of the loop —
not 10×.** The reason: at $n = 3$, `std::sort` takes only 0.23 branch misses per
push. **It was never mispredict-bound**, so removing the branches cannot buy what
that assumption predicted. Weighted by §1's shares this is ~7% of the campaign at
$10^7$ and ~5% at $10^9$, against a change that re-baselines the catalog's point
order. **Not recommended.**

### 2.4 Threads (B5) and the factor base (B6)

Thread scaling at $B = 10^8$, and the first measurement of the `alignas(64)`
false-sharing fix — which nothing could measure before, because ThreadSanitizer
cannot see false sharing and a single-threaded profile cannot either:

| threads | Mnorm/s | speedup | `alignas(64)` gain |
| --- | --- | --- | --- |
| 1 | 8.99 | 1.00× | **0.2%** |
| 2 | 17.82 | 1.98× | 4.6% |
| 4 | 35.03 | 3.89× | 4.0% |
| 8 | 68.76 | 7.64× | 4.2% |
| 16 | 84.37 | 9.38× | 3.4% |

Nothing at one thread, ~4% at two or more: **that is exactly the false-sharing
signature**, and it confirms the fix was worth making. Scaling is 7.64× on 8
cores; SMT adds only 1.23× beyond that, as expected for ALU-bound work.

**The factor base is not a cache problem** (B6, $B = 10^{10}$, 33 MB of records
against a 16 MB L3):

| | cyc/norm | L1d misses/norm |
| --- | --- | --- |
| `Map` (mmap) | 245.3 | 0.0045 |
| `Load` (heap) | 244.8 | 0.0409 |

**`Map` costs +0.2% cycles against `Load` in steady state** — this settles §9.4's
deliberately-unanswered question, because cycles are frequency-invariant where
its wall clock was not. And at 0.005–0.04 L1 misses per norm the walk essentially
never misses: the DFS strides the base ascending, and the hardware prefetcher
handles it. The 2–3× cache pessimism budgeted in §3 for $10^{12}$ looks too
pessimistic.

---

## 3. GPU primitive ceilings

All measured with clock discipline; the achieved clock is in every row of the
raw output. Nominal INT32 roofline is 24 SM × 64 lanes × clock.

| | primitive | measured | note |
| --- | --- | --- | --- |
| **P1** | `uint32` add | 5.13e12 op/s | 177% of nominal — simple adds evidently issue beyond the INT32-designated lanes on Ada |
| | `uint32` multiply | 2.89e12 op/s | **99.4% of nominal** — IMAD really is 64/SM/clk, which validates the model |
| | `uint64` add | 2.72e12 op/s | counted as 2×32-bit above |
| | **meet inner loop** | **7.34e11 quad-partitions/s** | 3 partitions from registers, the number a device meet must beat |
| **P2** | atomic append, hit rate 0 | 1.69e11 thread/s | |
| | atomic append, hit rate 1e-6 | **1.69e11 thread/s** | *identical* — the campaign's operating point (88 hits in 165 M norms) is free |
| | atomic append, hit rate 1 | 4.18e10 thread/s | 4× slower only when every thread appends |
| **P3** | CUB radix sort, 1.34e8 pairs | 8.4e8 key/s, 54 GB/s | **S-E's 1.374e11 pairs: 163 s** |
| **P4** | `atomicCAS` hash insert | 1.0e9 insert/s | **S-E's 1.374e11 inserts: 138 s**, including the table clear |
| **P5** | bitset rotate (funnelshift) | 1.5e10 → 1.1e11 word-op/s | $p$ = 1024 → 65536; small masks do not fill the machine |
| **P6** | 256-bit compare + add | 1.02e11 op-pair/s | **S-E's exact windowing: 1.4 s** |
| **P7** | `cudaHostAlloc` per transfer | 1.45 GB/s | the current code path |
| | pooled pinned buffer | **13.4 GB/s** | **9.2× faster**, 84% of PCIe 4.0 ×8 |
| | pageable | 13.1 GB/s | the win is not pinning, it is not allocating |

---

## 4. What the ceilings decide

**S-E (milestone G4) is confirmed, and cheaper than the plan projected.** The
plan estimated "minutes"; the measurement says the collision detection is
**138 s by hash insert, 163 s by radix sort**, and the exact 256-bit windowing
that §3.3 flagged as an uncosted risk is **1.4 s — free**. Hashing narrowly beats
sorting and needs a quarter of the memory, so build P4's shape, not P3's. The
CPU alternative pays a cache miss per key on a table far larger than L3; the GPU
advantage here is real and this is the search to build.

**S-D's transition is GPU-shaped but not yet worth building.** P5 confirms the
funnelshift rotate runs at 1.1e11 word-ops/s at large $p$ — but only 1.5e10 at
$p = 1024$, where the mask is too small to fill the device. Since the design
target is $p \le 10^3$, **the GPU is at its worst exactly where S-D lives.** That
strengthens rev. 3's "CPU first" verdict on independent grounds, and it means the
GPU version is worth writing only if the local-model question (§2.2 of the plan,
still the top open item) pushes $p$ past $10^4$.

**The `cudaHostAlloc` defect is a measured 9.2× on H2D**, not a suspicion. Any
sized run needs the pooled session object first.

**The atomic-append design assumption holds.** Design §4.3's "the atomic is
uncontended and the write bandwidth is unmeasurable" is confirmed: hit rate 1e-6
is *indistinguishable* from hit rate 0. Do not design for it.

**S-A on device: the arithmetic was never the problem.** The GPU meet inner loop
runs at 7.34e11 quad-partitions/s against roughly 1.6e10 for 16 CPU threads
(0.36–0.66 ns/quad-partition measured single-thread, scaled by B5's 9.38×) — a
**~46× ceiling**. But the design bound's *bucketed* meet work is 1.48e12
partitions, i.e. **~2 seconds of device arithmetic**. The meet is not what makes
S-A expensive; generation is. This is the same conclusion rev. 3 reached from
cost models, now with a number: porting the meet alone would buy nothing.

---

## 5. What is still unmeasurable here

- **`ncu` is blocked** — `ERR_NVGPUCTRPERM`, driver 550 defaulting to
  `NVreg_RestrictProfilingToAdminUsers=1`. Warp stall reasons and memory-pipe
  utilisation are unavailable, so §3 says what the device achieves but not why it
  falls short. To enable:
  ```sh
  echo 'options nvidia NVreg_RestrictProfilingToAdminUsers=0' |
    sudo tee /etc/modprobe.d/nvidia-profiling.conf
  sudo update-initramfs -u && reboot     # or just run: sudo ncu ...
  ```
- **The GPU clock sags under load** — 1.89 GHz on ALU work, 1.09–1.37 GHz during
  large sorts. P3's figures are therefore a *laptop* result; a desktop 4060 would
  do better and the S-E projection is conservative.
- **B5 reconstructs the campaign's `ThreadState`** rather than instantiating the
  shipped one (88 B here against the real object's larger footprint). The
  false-sharing signature is unambiguous, but the 4% is for the reconstruction.
- **Multi-threaded phase shares are unmeasured.** §1 is a single-thread profile;
  how the phases scale relative to each other across 16 threads is not known.
