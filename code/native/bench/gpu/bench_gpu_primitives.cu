// P1-P7 -- hardware ceilings for the primitives the GPU milestones need.
//
// These are not kernels for the searches. They are the question "what can this
// device do at all?", asked once per primitive, so that
// drafts/gpu_implementation_plan.md's cost models can be replaced with
// measurements. Read gpu_harness.cuh first: every number below is meaningless
// without the clock printed beside it.
//
//   P1  integer ALU ceiling            is "the arithmetic is free" true?
//   P2  atomic append throughput       DeviceAppendBuffer, at real hit rates
//   P3  CUB radix sort                 S-E's critical path
//   P4  atomicCAS hash insert          S-E's alternative, and S-D's dedup
//   P5  bitset rotate (funnelshift)    S-D's transition
//   P6  256-bit compare and add        S-E's exact windowing
//   P7  pinned staging and H2D         the known-unfixed pooling defect
//   P8  the same add, three ways       limb width and the inline PTX carry chain
//
//   ./bench_gpu_primitives [p1 p2 ...]      (default: all)
#include "bench/gpu/gpu_harness.cuh"

#include "core/uint.h"

#include <cub/device/device_radix_sort.cuh>

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

using namespace pte::gpubench;

namespace {

constexpr int kBlock = 256;

// S-E enumerates every unordered pair of a champion norm's D = 2^19 points:
// C(524288, 2) = 1.374e11. Value-windowing enumerates each pair exactly once,
// so this total is what every S-E projection below divides by a measured rate.
constexpr double kSEPairs = 1.374e11;

// ===========================================================================
// P1 -- integer ALU ceiling
// ===========================================================================
//
// Eight independent accumulators, so the loop is throughput-bound rather than
// latency-bound. A dependent chain would measure ALU *latency* and report a
// number 4x too low, which is one classic way to "prove" a GPU is slow at
// integer work.
//
// The other classic way is to get it 12x too *high*, which is what a first
// draft of this file did: written as plain `a += constant` in a loop, nvcc
// strength-reduces the whole thing to `a += constant * reps` and the benchmark
// reports 1231% of the hardware roofline. Any ALU measurement exceeding 100% of
// roofline is not a discovery, it is a deleted loop.
//
// So the operations are emitted as `asm volatile` PTX. The compiler may still
// schedule them, but it may not remove them or fold them together.

__global__ void AluAdd32(std::uint32_t *out, int reps) {
  std::uint32_t a = threadIdx.x, b = a + 1, c = a + 2, d = a + 3;
  std::uint32_t e = a + 4, f = a + 5, g = a + 6, h = a + 7;
  const std::uint32_t k = blockIdx.x | 1u;
  for (int i = 0; i < reps; ++i) {
    asm volatile("add.u32 %0, %0, %1;" : "+r"(a) : "r"(k));
    asm volatile("add.u32 %0, %0, %1;" : "+r"(b) : "r"(k));
    asm volatile("add.u32 %0, %0, %1;" : "+r"(c) : "r"(k));
    asm volatile("add.u32 %0, %0, %1;" : "+r"(d) : "r"(k));
    asm volatile("add.u32 %0, %0, %1;" : "+r"(e) : "r"(k));
    asm volatile("add.u32 %0, %0, %1;" : "+r"(f) : "r"(k));
    asm volatile("add.u32 %0, %0, %1;" : "+r"(g) : "r"(k));
    asm volatile("add.u32 %0, %0, %1;" : "+r"(h) : "r"(k));
  }
  if (threadIdx.x == 0xffff)
    out[blockIdx.x] = a + b + c + d + e + f + g + h;
}

__global__ void AluAdd64(std::uint64_t *out, int reps) {
  std::uint64_t a = threadIdx.x, b = a + 1, c = a + 2, d = a + 3;
  std::uint64_t e = a + 4, f = a + 5, g = a + 6, h = a + 7;
  const std::uint64_t k = blockIdx.x | 1u;
  for (int i = 0; i < reps; ++i) {
    asm volatile("add.u64 %0, %0, %1;" : "+l"(a) : "l"(k));
    asm volatile("add.u64 %0, %0, %1;" : "+l"(b) : "l"(k));
    asm volatile("add.u64 %0, %0, %1;" : "+l"(c) : "l"(k));
    asm volatile("add.u64 %0, %0, %1;" : "+l"(d) : "l"(k));
    asm volatile("add.u64 %0, %0, %1;" : "+l"(e) : "l"(k));
    asm volatile("add.u64 %0, %0, %1;" : "+l"(f) : "l"(k));
    asm volatile("add.u64 %0, %0, %1;" : "+l"(g) : "l"(k));
    asm volatile("add.u64 %0, %0, %1;" : "+l"(h) : "l"(k));
  }
  if (threadIdx.x == 0xffff)
    out[blockIdx.x] = a + b + c + d + e + f + g + h;
}

// 32-bit multiply, the hoist's operation (nu * nu).
__global__ void AluMul32(std::uint32_t *out, int reps) {
  std::uint32_t a = threadIdx.x | 1u, b = a + 2, c = a + 4, d = a + 6;
  const std::uint32_t k = (blockIdx.x | 1u) * 2654435761u;
  for (int i = 0; i < reps; ++i) {
    asm volatile("mul.lo.u32 %0, %0, %1;" : "+r"(a) : "r"(k));
    asm volatile("mul.lo.u32 %0, %0, %1;" : "+r"(b) : "r"(k));
    asm volatile("mul.lo.u32 %0, %0, %1;" : "+r"(c) : "r"(k));
    asm volatile("mul.lo.u32 %0, %0, %1;" : "+r"(d) : "r"(k));
  }
  if (threadIdx.x == 0xffff)
    out[blockIdx.x] = a + b + c + d;
}

// The meet's actual inner operation: sq[a] + sq[b] == sq[c] + sq[d], with the
// squares already hoisted into registers. Two adds and a compare on uint64.
__global__ void MeetKernelShape(std::uint32_t *out, int reps) {
  std::uint64_t s0 = threadIdx.x * 2654435761ull + 1;
  std::uint64_t s1 = s0 * 3, s2 = s0 * 5, s3 = s0 * 7;
  std::uint32_t hits = 0;
  for (int i = 0; i < reps; ++i) {
    hits += (s0 + s1 == s2 + s3) ? 1u : 0u;
    hits += (s0 + s2 == s1 + s3) ? 1u : 0u;
    hits += (s0 + s3 == s1 + s2) ? 1u : 0u;
    s0 += 1; s1 += 3; s2 += 5; s3 += 7;
  }
  if (threadIdx.x == 0xffff)
    out[blockIdx.x] = hits;
}

void P1(const Device &dev) {
  const int blocks = dev.prop.multiProcessorCount * 16;
  const int reps = 4096;
  const double threads = static_cast<double>(blocks) * kBlock;
  std::uint32_t *d32 = nullptr;
  std::uint64_t *d64 = nullptr;
  PTE_CUDA_OK(cudaMalloc(&d32, blocks * sizeof(std::uint32_t)));
  PTE_CUDA_OK(cudaMalloc(&d64, blocks * sizeof(std::uint64_t)));

  Table t("P1 integer ALU ceiling", "op/s");
  {
    const Measurement m =
        Time(dev, 20, [&] { AluAdd32<<<blocks, kBlock>>>(d32, reps); });
    t.Row("uint32 add, 8 independent chains", m, threads * reps * 8,
          dev.Int32Roofline(m.ghz));
  }
  {
    const Measurement m =
        Time(dev, 20, [&] { AluAdd64<<<blocks, kBlock>>>(d64, reps); });
    // A 64-bit add is two 32-bit ops (add + addc), so count it as two.
    t.Row("uint64 add (counted as 2 x 32-bit)", m, threads * reps * 8 * 2,
          dev.Int32Roofline(m.ghz));
  }
  {
    const Measurement m =
        Time(dev, 20, [&] { AluMul32<<<blocks, kBlock>>>(d32, reps); });
    t.Row("uint32 multiply (the sq[] hoist)", m, threads * reps * 4,
          dev.Int32Roofline(m.ghz));
  }
  {
    const Measurement m =
        Time(dev, 20, [&] { MeetKernelShape<<<blocks, kBlock>>>(d32, reps); });
    t.Row("meet inner loop, 3 partitions", m, threads * reps * 3, 0);
    t.Note("^ that rate is quad-partitions/s: the number to beat when a device "
           "meet is written");
  }
  PTE_CUDA_OK(cudaFree(d32));
  PTE_CUDA_OK(cudaFree(d64));
}

// ===========================================================================
// P2 -- atomic append at realistic hit rates
// ===========================================================================

struct AppendView {
  std::uint64_t *data;
  unsigned long long *counter;
  unsigned long long capacity;
};

__global__ void AppendKernel(AppendView out, std::uint32_t threshold, int n) {
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= n)
    return;
  // splitmix32, so "does this thread append" is deterministic and uncorrelated
  // with lane index -- a threshold on tid alone would make whole warps append
  // together and measure a different contention pattern than the real one.
  std::uint32_t h = static_cast<std::uint32_t>(tid) * 0x9e3779b9u;
  h ^= h >> 16;
  h *= 0x85ebca6bu;
  h ^= h >> 13;
  if (h < threshold) {
    const unsigned long long slot = atomicAdd(out.counter, 1ull);
    if (slot < out.capacity)
      out.data[slot] = tid;
  }
}

void P2(const Device &dev) {
  const int n = 1 << 24; // 16.8 M threads, ~ a campaign batch
  const int blocks = (n + kBlock - 1) / kBlock;
  const std::uint64_t cap = 1 << 22;
  std::uint64_t *data = nullptr;
  unsigned long long *counter = nullptr;
  PTE_CUDA_OK(cudaMalloc(&data, cap * sizeof(std::uint64_t)));
  PTE_CUDA_OK(cudaMalloc(&counter, sizeof(unsigned long long)));

  Table t("P2 atomic append (DeviceAppendBuffer::push shape)", "thread/s");
  const double rates[] = {0.0, 1e-6, 1e-3, 1e-2, 1e-1, 1.0};
  for (const double rate : rates) {
    const auto threshold =
        static_cast<std::uint32_t>(rate * 4294967295.0);
    PTE_CUDA_OK(cudaMemset(counter, 0, sizeof(unsigned long long)));
    AppendView v{data, counter, cap};
    const Measurement m = Time(dev, 20, [&] {
      AppendKernel<<<blocks, kBlock>>>(v, threshold, n);
    });
    char name[80];
    std::snprintf(name, sizeof(name), "hit rate %-9.6g%s", rate,
                  rate == 1e-6 ? "  <- the real campaign" : "");
    t.Row(name, m, n, 0);
  }
  t.Note("the campaign found 88 candidates in 165 M norms = 5e-7, so the top "
         "rows are the operating point");
  PTE_CUDA_OK(cudaFree(data));
  PTE_CUDA_OK(cudaFree(counter));
}

// ===========================================================================
// P3 -- CUB radix sort: S-E's critical path
// ===========================================================================

__global__ void FillKeys(std::uint64_t *k, std::uint64_t *v, int n,
                         std::uint64_t seed) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n)
    return;
  std::uint64_t x = seed + static_cast<std::uint64_t>(i);
  x ^= x >> 33; x *= 0xff51afd7ed558ccdull;
  x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ull;
  x ^= x >> 33;
  k[i] = x;
  v[i] = static_cast<std::uint64_t>(i);
}

void P3(const Device &dev) {
  std::size_t free_bytes = 0, total_bytes = 0;
  PTE_CUDA_OK(cudaMemGetInfo(&free_bytes, &total_bytes));

  Table t("P3 CUB radix sort, 8-byte key + 8-byte value", "key/s");
  t.Note("S-E windows hold ~1.4e8 pairs each; " +
         std::to_string(free_bytes / 1048576) + " MB free on device");

  for (const int n : {1 << 20, 1 << 23, 1 << 25, 1 << 26, 1 << 27}) {
    // in + out for keys and values, plus CUB temp storage.
    const std::size_t need =
        static_cast<std::size_t>(n) * 8 * 4 + (static_cast<std::size_t>(n) * 8);
    if (need > free_bytes * 8 / 10) {
      char name[80];
      std::snprintf(name, sizeof(name), "n = %-12d skipped, needs %zu MB", n,
                    need / 1048576);
      t.Note(name);
      continue;
    }
    std::uint64_t *ki = nullptr, *ko = nullptr, *vi = nullptr, *vo = nullptr;
    PTE_CUDA_OK(cudaMalloc(&ki, n * sizeof(std::uint64_t)));
    PTE_CUDA_OK(cudaMalloc(&ko, n * sizeof(std::uint64_t)));
    PTE_CUDA_OK(cudaMalloc(&vi, n * sizeof(std::uint64_t)));
    PTE_CUDA_OK(cudaMalloc(&vo, n * sizeof(std::uint64_t)));
    FillKeys<<<(n + kBlock - 1) / kBlock, kBlock>>>(ki, vi, n, 12345);
    PTE_CUDA_OK(cudaDeviceSynchronize());

    void *temp = nullptr;
    std::size_t temp_bytes = 0;
    cub::DeviceRadixSort::SortPairs(temp, temp_bytes, ki, ko, vi, vo, n);
    PTE_CUDA_OK(cudaMalloc(&temp, temp_bytes));

    const Measurement m = Time(dev, 8, [&] {
      cub::DeviceRadixSort::SortPairs(temp, temp_bytes, ki, ko, vi, vo, n);
    });
    char name[80];
    std::snprintf(name, sizeof(name), "n = %-12d (%zu MB temp)", n,
                  temp_bytes / 1048576);
    // A radix sort over 64-bit keys is several full read+write passes; report
    // against memory bandwidth, which is what actually binds it.
    t.Row(name, m, n, 0);
    // S-E's total is C(2^19, 2) = 1.374e11 pairs, and every pair is sorted
    // exactly once no matter how the windows are cut (that is the point of
    // windowing by value -- see gpu_implementation_plan.md 3.1). So the
    // projection is total/rate, NOT per-window-cost x window-count; a first
    // draft of this line multiplied by 4000 windows and reported 4x the truth.
    const double rate = static_cast<double>(n) / m.seconds;
    std::printf("      -> %.1f GB/s effective (16 B/key, 2 passes assumed); "
                "S-E's 1.374e11 pairs: %.0f s\n",
                static_cast<double>(n) * 16 * 2 * 2 / m.seconds / 1e9,
                kSEPairs / rate);

    PTE_CUDA_OK(cudaFree(ki)); PTE_CUDA_OK(cudaFree(ko));
    PTE_CUDA_OK(cudaFree(vi)); PTE_CUDA_OK(cudaFree(vo));
    PTE_CUDA_OK(cudaFree(temp));
  }
}

// ===========================================================================
// P4 -- open-addressing hash insert with atomicCAS
// ===========================================================================

__global__ void HashInsert(std::uint64_t *table, std::uint64_t mask, int n,
                           std::uint64_t seed, unsigned long long *collisions) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n)
    return;
  std::uint64_t x = seed + static_cast<std::uint64_t>(i);
  x ^= x >> 33; x *= 0xff51afd7ed558ccdull;
  x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ull;
  x ^= x >> 33;
  const std::uint64_t key = x | 1ull; // 0 marks an empty slot
  std::uint64_t slot = x & mask;
  for (int probe = 0; probe < 32; ++probe) {
    const std::uint64_t prev = atomicCAS(
        reinterpret_cast<unsigned long long *>(table + slot), 0ull,
        static_cast<unsigned long long>(key));
    if (prev == 0)
      return; // inserted
    if (prev == key) {
      atomicAdd(collisions, 1ull); // an equal key: exactly what S-E looks for
      return;
    }
    slot = (slot + 1) & mask;
  }
}

void P4(const Device &dev) {
  Table t("P4 atomicCAS open-addressing insert", "insert/s");
  const int table_log = 26; // 64 Mi slots = 512 MB
  const std::uint64_t slots = 1ull << table_log;
  std::uint64_t *table = nullptr;
  unsigned long long *coll = nullptr;
  PTE_CUDA_OK(cudaMalloc(&table, slots * sizeof(std::uint64_t)));
  PTE_CUDA_OK(cudaMalloc(&coll, sizeof(unsigned long long)));

  for (const double load : {0.1, 0.25, 0.5, 0.75}) {
    const int n = static_cast<int>(static_cast<double>(slots) * load);
    const int blocks = (n + kBlock - 1) / kBlock;
    const Measurement m = Time(dev, 8, [&] {
      cudaMemsetAsync(table, 0, slots * sizeof(std::uint64_t));
      cudaMemsetAsync(coll, 0, sizeof(unsigned long long));
      HashInsert<<<blocks, kBlock>>>(table, slots - 1, n, 999, coll);
    });
    char name[80];
    std::snprintf(name, sizeof(name), "load factor %.2f, n = %d", load, n);
    t.Row(name, m, n, 0);
    std::printf("      -> S-E's 1.374e11 inserts: %.0f s\n",
                kSEPairs * m.seconds / n);
  }
  t.Note("includes the table clear, which is what a per-window S-E pass pays "
         "too");
  PTE_CUDA_OK(cudaFree(table));
  PTE_CUDA_OK(cudaFree(coll));
}

// ===========================================================================
// P5 -- bitset rotate: S-D's transition
// ===========================================================================
//
// T' = { t^2 - L : t in A } is, once the squares of A are known as a mask S, a
// cyclic shift of S by L. __funnelshift_l gives the cross-word shift in one
// instruction; this measures how fast a p-bit mask can be rotated.

__global__ void RotateMask(const std::uint32_t *__restrict__ src,
                           std::uint32_t *__restrict__ dst, int words,
                           int shifts) {
  extern __shared__ std::uint32_t sh[];
  for (int i = threadIdx.x; i < words; i += blockDim.x)
    sh[i] = src[i];
  __syncthreads();

  std::uint32_t acc = 0;
  for (int s = 0; s < shifts; ++s) {
    const int wordshift = s / 32, bitshift = s % 32;
    for (int i = threadIdx.x; i < words; i += blockDim.x) {
      const int lo = (i - wordshift + words) % words;
      const int hi = (lo - 1 + words) % words;
      acc ^= __funnelshift_l(sh[hi], sh[lo], bitshift);
    }
  }
  if (threadIdx.x == 0xffff)
    dst[blockIdx.x] = acc;
}

void P5(const Device &dev) {
  Table t("P5 bitset rotate (S-D transition)", "word-op/s");
  for (const int p : {1024, 8192, 65536}) {
    const int words = p / 32;
    const int shifts = 256;
    const int blocks = dev.prop.multiProcessorCount * 8;
    std::uint32_t *src = nullptr, *dst = nullptr;
    PTE_CUDA_OK(cudaMalloc(&src, words * sizeof(std::uint32_t)));
    PTE_CUDA_OK(cudaMalloc(&dst, blocks * sizeof(std::uint32_t)));
    PTE_CUDA_OK(cudaMemset(src, 0x5a, words * sizeof(std::uint32_t)));
    const std::size_t shmem = words * sizeof(std::uint32_t);
    if (shmem > dev.prop.sharedMemPerBlock) {
      char name[80];
      std::snprintf(name, sizeof(name), "p = %-8d skipped: %zu B > shared %zu B",
                    p, shmem, dev.prop.sharedMemPerBlock);
      t.Note(name);
      PTE_CUDA_OK(cudaFree(src));
      PTE_CUDA_OK(cudaFree(dst));
      continue;
    }
    const Measurement m = Time(dev, 15, [&] {
      RotateMask<<<blocks, kBlock, shmem>>>(src, dst, words, shifts);
    });
    char name[80];
    std::snprintf(name, sizeof(name), "p = %-6d (%d words, %d shifts/block)", p,
                  words, shifts);
    t.Row(name, m, static_cast<double>(blocks) * words * shifts, 0);
    PTE_CUDA_OK(cudaFree(src));
    PTE_CUDA_OK(cudaFree(dst));
  }
  t.Note("S-D at p = 1000 needs ~p^4/64 = 1.5e10 word-ops per prime");
}

// ===========================================================================
// P6 -- 256-bit compare and add: S-E's exact windowing
// ===========================================================================

struct U256 {
  std::uint64_t w[4]; // little-endian limbs
};

__device__ inline int Cmp256(const U256 &a, const U256 &b) {
  for (int i = 3; i >= 0; --i) {
    if (a.w[i] != b.w[i])
      return a.w[i] < b.w[i] ? -1 : 1;
  }
  return 0;
}

__device__ inline U256 Add256(const U256 &a, const U256 &b) {
  U256 r;
  unsigned long long carry = 0;
  for (int i = 0; i < 4; ++i) {
    const unsigned long long s = a.w[i] + b.w[i];
    const unsigned long long c1 = s < a.w[i] ? 1ull : 0ull;
    const unsigned long long s2 = s + carry;
    carry = c1 | (s2 < s ? 1ull : 0ull);
    r.w[i] = s2;
  }
  return r;
}

__global__ void Bench256(std::uint32_t *out, int reps) {
  U256 a, b;
  for (int i = 0; i < 4; ++i) {
    a.w[i] = threadIdx.x * 0x9e3779b97f4a7c15ull + i;
    b.w[i] = threadIdx.x * 0xc2b2ae3512345679ull + i;
  }
  std::uint32_t acc = 0;
  for (int i = 0; i < reps; ++i) {
    a = Add256(a, b);
    acc += static_cast<std::uint32_t>(Cmp256(a, b) + 1);
  }
  if (threadIdx.x == 0xffff)
    out[blockIdx.x] = acc + static_cast<std::uint32_t>(a.w[0]);
}

void P6(const Device &dev) {
  const int blocks = dev.prop.multiProcessorCount * 16;
  const int reps = 2048;
  std::uint32_t *out = nullptr;
  PTE_CUDA_OK(cudaMalloc(&out, blocks * sizeof(std::uint32_t)));
  Table t("P6 256-bit compare + add (S-E windowing)", "op-pair/s");
  const Measurement m =
      Time(dev, 20, [&] { Bench256<<<blocks, kBlock>>>(out, reps); });
  const double ops = static_cast<double>(blocks) * kBlock * reps;
  t.Row("4-limb add + lexicographic compare", m, ops, 0);
  t.Note("S-E enumerates 1.4e11 pairs; at this rate the exact windowing "
         "comparison alone costs the time shown x 1.4e11/ops");
  std::printf("      -> %.1f s for S-E's 1.4e11 exact comparisons\n",
              1.4e11 / (ops / m.seconds));
  PTE_CUDA_OK(cudaFree(out));
}

// ===========================================================================
// P7 -- pinned staging: the known-unfixed pooling defect
// ===========================================================================

void P7(const Device &dev) {
  Table t("P7 pinned host staging and H2D", "byte/s");
  const std::size_t bytes = 64u << 20; // 64 MB, a plausible factor-base segment
  void *device = nullptr;
  PTE_CUDA_OK(cudaMalloc(&device, bytes));

  // (a) allocate a fresh pinned buffer every transfer -- what the demo does.
  {
    const Measurement m = Time(dev, 5, [&] {
      void *host = nullptr;
      cudaHostAlloc(&host, bytes, cudaHostAllocDefault);
      cudaMemcpy(device, host, bytes, cudaMemcpyHostToDevice);
      cudaFreeHost(host);
    });
    t.Row("cudaHostAlloc per transfer (current)", m, bytes,
          dev.MemoryRoofline());
  }
  // (b) allocate once, reuse -- what a session object would do.
  {
    void *host = nullptr;
    PTE_CUDA_OK(cudaHostAlloc(&host, bytes, cudaHostAllocDefault));
    std::memset(host, 0, bytes);
    const Measurement m = Time(dev, 20, [&] {
      cudaMemcpy(device, host, bytes, cudaMemcpyHostToDevice);
    });
    t.Row("pooled pinned buffer (proposed)", m, bytes, 0);
    PTE_CUDA_OK(cudaFreeHost(host));
  }
  // (c) pageable, for reference.
  {
    std::vector<char> host(bytes, 0);
    const Measurement m = Time(dev, 10, [&] {
      cudaMemcpy(device, host.data(), bytes, cudaMemcpyHostToDevice);
    });
    t.Row("pageable host memory", m, bytes, 0);
  }
  t.Note("PCIe 4.0 x8 is ~16 GB/s theoretical; the whole 2.7 GB factor base "
         "crosses in well under a second either way");
  PTE_CUDA_OK(cudaFree(device));
}

// ===========================================================================
// P8 -- the same 256-bit add + compare, three ways
// ===========================================================================
//
// P6 answered "can we afford exact 256-bit windowing?" (yes, 1.4 s of S-E's
// 138). P8 answers the follow-up that decides how core/uint.h is written:
// **does the limb width or the carry mechanism matter, and by how much?**
//
// Three implementations of one operation:
//
//   (a) 4 x uint64 limbs, portable C   -- P6's own U256, re-run here so the two
//                                         tables are comparable and so a
//                                         disagreement would be visible
//   (b) 8 x uint32 limbs, portable C   -- the documented limb convention,
//                                         carry rebuilt by 64-bit accumulate
//   (c) 8 x uint32 limbs, inline PTX   -- add.cc.u32 / addc.cc.u32 / addc.u32,
//                                         i.e. core/uint.h as it ships
//
// The loop shape is P6's: one dependent chain, so this is *latency*, which is
// the honest number for a carry chain (a carry chain is a dependency by
// definition -- there is no independent-accumulator version of it). Row (a)
// reproducing P6 is the check that the harness is measuring the same thing.
//
// Verify the SASS alongside the timing, since that is the actual claim:
//   cuobjdump -sass build/dev/bench/gpu/bench_gpu_primitives |
//     grep -A40 Bench256U32Ptx

struct U256x32 {
  std::uint32_t w[8];
};

__device__ inline int Cmp256x32(const U256x32 &a, const U256x32 &b) {
  for (int i = 7; i >= 0; --i) {
    if (a.w[i] != b.w[i])
      return a.w[i] < b.w[i] ? -1 : 1;
  }
  return 0;
}

// (b) the portable form: accumulate in 64 bits and shift the carry out.
__device__ inline U256x32 Add256x32C(const U256x32 &a, const U256x32 &b) {
  U256x32 r;
  std::uint64_t carry = 0;
#pragma unroll
  for (int i = 0; i < 8; ++i) {
    const std::uint64_t s = static_cast<std::uint64_t>(a.w[i]) + b.w[i] + carry;
    r.w[i] = static_cast<std::uint32_t>(s);
    carry = s >> 32;
  }
  return r;
}

__global__ void Bench256U64(std::uint32_t *out, int reps) {
  U256 a, b;
  for (int i = 0; i < 4; ++i) {
    a.w[i] = threadIdx.x * 0x9e3779b97f4a7c15ull + i;
    b.w[i] = threadIdx.x * 0xc2b2ae3512345679ull + i;
  }
  std::uint32_t acc = 0;
  for (int i = 0; i < reps; ++i) {
    a = Add256(a, b);
    acc += static_cast<std::uint32_t>(Cmp256(a, b) + 1);
  }
  if (threadIdx.x == 0xffff)
    out[blockIdx.x] = acc + static_cast<std::uint32_t>(a.w[0]);
}

__global__ void Bench256U32C(std::uint32_t *out, int reps) {
  U256x32 a, b;
  for (int i = 0; i < 8; ++i) {
    a.w[i] = threadIdx.x * 0x9e3779b9u + i;
    b.w[i] = threadIdx.x * 0x12345679u + i;
  }
  std::uint32_t acc = 0;
  for (int i = 0; i < reps; ++i) {
    a = Add256x32C(a, b);
    acc += static_cast<std::uint32_t>(Cmp256x32(a, b) + 1);
  }
  if (threadIdx.x == 0xffff)
    out[blockIdx.x] = acc + a.w[0];
}

// (c) core/uint.h itself -- the shipping implementation, not a copy of it, so
// this row cannot drift away from the code it is quoted at.
__global__ void Bench256U32Ptx(std::uint32_t *out, int reps) {
  pte::core::U256 a, b;
  for (int i = 0; i < 8; ++i) {
    a.w[i] = threadIdx.x * 0x9e3779b9u + i;
    b.w[i] = threadIdx.x * 0x12345679u + i;
  }
  std::uint32_t acc = 0;
  for (int i = 0; i < reps; ++i) {
    a = pte::core::add(a, b);
    acc += static_cast<std::uint32_t>(pte::core::cmp(a, b) + 1);
  }
  if (threadIdx.x == 0xffff)
    out[blockIdx.x] = acc + a.w[0];
}

void P8(const Device &dev) {
  const int blocks = dev.prop.multiProcessorCount * 16;
  const int reps = 2048;
  std::uint32_t *out = nullptr;
  PTE_CUDA_OK(cudaMalloc(&out, blocks * sizeof(std::uint32_t)));
  const double ops = static_cast<double>(blocks) * kBlock * reps;

  Table t("P8 256-bit add + compare, by limb width and carry", "op-pair/s");
  const auto row = [&](const char *name, const Measurement &m) {
    t.Row(name, m, ops, 0);
    std::printf("      -> %.2f s for S-E's %.3e add+compare pairs\n",
                kSEPairs / (ops / m.seconds), kSEPairs);
  };
  row("4 x uint64 limbs, portable C  (= P6)",
      Time(dev, 20, [&] { Bench256U64<<<blocks, kBlock>>>(out, reps); }));
  row("8 x uint32 limbs, portable C",
      Time(dev, 20, [&] { Bench256U32C<<<blocks, kBlock>>>(out, reps); }));
  row("8 x uint32 limbs, inline PTX  (core/uint.h)",
      Time(dev, 20, [&] { Bench256U32Ptx<<<blocks, kBlock>>>(out, reps); }));
  t.Note("one dependent chain, so this is carry-chain latency -- which is the "
         "only thing a carry chain has");
  PTE_CUDA_OK(cudaFree(out));
}

} // namespace

int main(int argc, char **argv) {
  std::vector<std::string> want(argv + 1, argv + argc);
  const auto run = [&](const char *tag) {
    return want.empty() ||
           std::find(want.begin(), want.end(), tag) != want.end();
  };

  const Device dev = Device::Query();
  dev.Print();
  std::printf("idle SM clock right now: %.2f GHz (max %.2f) -- every table "
              "below prints the clock it achieved\n",
              dev.CurrentGHz(), dev.MaxClockGHz());

  if (run("p1")) P1(dev);
  if (run("p2")) P2(dev);
  if (run("p3")) P3(dev);
  if (run("p4")) P4(dev);
  if (run("p5")) P5(dev);
  if (run("p6")) P6(dev);
  if (run("p7")) P7(dev);
  if (run("p8")) P8(dev);
  return 0;
}
