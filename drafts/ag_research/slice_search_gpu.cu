// C1 on the GPU: rational points of Bremner's slice to large height.
//
// Same search as slice_search.cpp (primitive (x:y:z), z >= 1, x >= y, |x|,|y|,z <= H,
// q1 = E + O and q2 = E - O both squares, degenerate lines/conics skipped), but the
// filter runs on the GPU: for each z the host builds, for 16 small primes p, a
// table over residues (x mod p, y mod p) marking where q1 and q2 are both squares
// or 0 mod p. Threads test their (x, y) against the tables (no wide arithmetic);
// survivors go back to the host, which does the exact __int128 square test.
// Build: nvcc -O3 -std=c++20 -ccbin=/usr/bin/g++-13 -arch=sm_89 -U_GNU_SOURCE -o slice_search_gpu slice_search_gpu.cu
//   (-U_GNU_SOURCE: glibc 2.43 declares cospi/sinpi noexcept and clashes with CUDA 12.9 headers)
// Run:   ./slice_search_gpu H [z_from z_to] > hits.txt
#include <cstdio>
#include <cstdlib>
#include <cstdint>

#include <numeric>
#include <vector>
#include <cuda_runtime.h>

typedef __int128 i128;
#define NP 16
#define ROWS 64
static const int PR[NP] = {5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61};

__constant__ int d_pr[NP];
__constant__ unsigned long long d_tab[NP * ROWS];   // row r of prime i: bitmask over y-residues
__constant__ int d_off[NP];                        // (H mod p): x + H, y + H shifted residues

#define CHECK(x) do { cudaError_t e = (x); if (e != cudaSuccess) { \
  fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(e), __FILE__, __LINE__); exit(1); } } while (0)

__global__ void sieve_kernel(long H, long z, unsigned long long cap, unsigned long long* count, int2* out) {
  long n = 2 * H + 1;
  unsigned long long idx = (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= (unsigned long long)n * n) return;
  long xs = idx / n, ys = idx % n;          // shifted coordinates in [0, 2H]
  long x = xs - H, y = ys - H;
  if (y > x) return;
  if (x == z || x == -z || y == z || y == -z || x + y == 0) return;
  for (int i = 0; i < NP; ++i) {
    int p = d_pr[i];
    int rx = (int)(xs % p), ry = (int)(ys % p);
    if (!((d_tab[i * ROWS + rx] >> ry) & 1ULL)) return;
  }
  long long c1 = y*y - x*y - x*z - y*z - 10*z*z, c2 = x*x - x*y - x*z - y*z - 10*z*z;
  long long c3 = x*x - x*y + x*z + y*z - 10*z*z, c4 = y*y - x*y + x*z + y*z - 10*z*z;
  if (c1 == 0 || c2 == 0 || c3 == 0 || c4 == 0) return;
  unsigned long long k = atomicAdd(count, 1ULL);
  if (k < cap) out[k] = make_int2((int)x, (int)y);
}

static void qmod(long x, long y, long z, int p, int& q1, int& q2) {
  long X = ((x % p) + p) % p, Y = ((y % p) + p) % p, Z = ((z % p) + p) % p;
  long X2 = X * X % p, Y2 = Y * Y % p, Z2 = Z * Z % p, Z3 = Z2 * Z % p, Z4 = Z2 * Z2 % p, Z5 = Z4 * Z % p;
  long E = (X2 * X2 % p * Y2 - 2 * X2 * X % p * Y2 % p * Y + X2 * Y2 % p * Y2
            + Z2 * ((10 * X2 * X % p * Y - 20 * X2 * Y2 + 10 * X * Y2 % p * Y) % p)
            + Z4 * ((45 * X2 - 10 * X * Y + 45 * Y2) % p)) % p;
  long O = (-4 * Z3 * ((X2 * X + X2 * Y + X * Y2 + Y2 * Y) % p) + 40 * Z5 * ((X + Y) % p)) % p;
  q1 = (int)(((E + O) % p + p) % p);
  q2 = (int)(((E - O) % p + p) % p);
}

static bool is_square(i128 v) {
  if (v < 0) return false;
  long double r = __builtin_sqrtl((long double)v);
  i128 s = (i128)r;
  for (i128 d = s - 2; d <= s + 2; ++d) if (d >= 0 && d * d == v) return true;
  return false;
}

static void q_exact(long x, long y, long z, i128& q1, i128& q2) {
  i128 X = x, Y = y, Z = z, X2 = X*X, X3 = X2*X, X4 = X2*X2, Y2 = Y*Y, Y3 = Y2*Y, Y4 = Y2*Y2;
  i128 Z2 = Z*Z, Z3 = Z2*Z, Z4 = Z2*Z2, Z5 = Z4*Z;
  i128 E = X4*Y2 - 2*X3*Y3 + X2*Y4 + Z2*(10*X3*Y - 20*X2*Y2 + 10*X*Y3) + Z4*(45*X2 - 10*X*Y + 45*Y2);
  i128 O = -4*Z3*(X3 + X2*Y + X*Y2 + Y3) + 40*Z5*(X + Y);
  q1 = E + O; q2 = E - O;
}

int main(int argc, char** argv) {
  long H = argc > 1 ? atol(argv[1]) : 1000;
  long z0 = argc > 2 ? atol(argv[2]) : 1, z1 = argc > 3 ? atol(argv[3]) : H;
  CHECK(cudaMemcpyToSymbol(d_pr, PR, sizeof(PR)));
  const unsigned long long cap = 1ULL << 24;
  unsigned long long* d_count; int2* d_out;
  CHECK(cudaMalloc(&d_count, sizeof(unsigned long long)));
  CHECK(cudaMalloc(&d_out, cap * sizeof(int2)));
  std::vector<int2> h_out(cap);
  long n = 2 * H + 1;
  unsigned long long total_threads = (unsigned long long)n * n;
  int block = 256;
  unsigned long long grid = (total_threads + block - 1) / block;
  long long nsurv = 0, nhits = 0;
  for (long z = z0; z <= z1; ++z) {
    // tables for this z, indexed by shifted residues (x + H) mod p, (y + H) mod p
    unsigned long long tab[NP * ROWS] = {0};
    for (int i = 0; i < NP; ++i) {
      int p = PR[i];
      std::vector<char> sq(p, 0);
      for (int a = 0; a < p; ++a) sq[(a * a) % p] = 1;
      for (int rx = 0; rx < p; ++rx)
        for (int ry = 0; ry < p; ++ry) {
          long x = rx - (H % p), y = ry - (H % p);     // any representatives
          int q1, q2; qmod(x, y, z, p, q1, q2);
          if (sq[q1] && sq[q2]) tab[i * ROWS + rx] |= (1ULL << ry);
        }
    }
    CHECK(cudaMemcpyToSymbol(d_tab, tab, sizeof(tab)));
    CHECK(cudaMemset(d_count, 0, sizeof(unsigned long long)));
    sieve_kernel<<<(unsigned)grid, block>>>(H, z, cap, d_count, d_out);
    CHECK(cudaGetLastError());
    unsigned long long cnt;
    CHECK(cudaMemcpy(&cnt, d_count, sizeof(cnt), cudaMemcpyDeviceToHost));
    if (cnt > cap) { fprintf(stderr, "survivor overflow at z = %ld (%llu)\n", z, cnt); exit(2); }
    CHECK(cudaMemcpy(h_out.data(), d_out, cnt * sizeof(int2), cudaMemcpyDeviceToHost));
    nsurv += cnt;
    for (unsigned long long k = 0; k < cnt; ++k) {
      long x = h_out[k].x, y = h_out[k].y;
      if (std::gcd(std::gcd(std::labs(x), std::labs(y)), z) != 1) continue;
      i128 q1, q2; q_exact(x, y, z, q1, q2);
      if (is_square(q1) && is_square(q2)) { printf("%ld %ld %ld\n", x, y, z); ++nhits; }
    }
    if (z % 500 == 0) { fprintf(stderr, "z = %ld: survivors so far %lld, hits %lld\n", z, nsurv, nhits); fflush(stdout); }
  }
  fprintf(stderr, "H = %ld, z in [%ld, %ld]: GPU survivors %lld, exact hits %lld\n", H, z0, z1, nsurv, nhits);
  return 0;
}
