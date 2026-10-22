// C1: rational points of Bremner's slice to large height (exact, __int128).
//
// Finds primitive (x:y:z), z >= 1, x >= y, |x|, |y|, z <= H, with
//   q1 = E + O  and  q2 = E - O  both perfect squares,
// where q1(x,y,z) is Bremner's sextic (2-5), q2(x,y,z) = q1(x,y,-z), and E/O are
// its parts even/odd in z. |q| <= 240 H^6 < 2^127 for H <= 2.5e5.
// Skipped as degenerate: z = 0, x = +-z, y = +-z, x + y = 0 (split lines) and the
// four L4 = 0 conics. Output: one line "x y z" per hit (to stdout).
// Build: g++ -O3 -march=native -fopenmp -o slice_search slice_search.cpp
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <numeric>

typedef __int128 i128;

static bool sq64[64], sq63[63], sq65[65], sq11[11];

static inline bool maybe_square(i128 v) {
  if (v < 0) return false;
  if (!sq64[(uint64_t)v & 63]) return false;
  if (!sq63[(int)(v % 63)]) return false;
  if (!sq65[(int)(v % 65)]) return false;
  if (!sq11[(int)(v % 11)]) return false;
  return true;
}

static inline bool is_square(i128 v) {
  if (!maybe_square(v)) return false;
  long double r = sqrtl((long double)v);
  i128 s = (i128)r;
  for (i128 d = s - 2; d <= s + 2; ++d)
    if (d >= 0 && d * d == v) return true;
  return false;
}

static inline bool degenerate(long x, long y, long z) {
  if (x == z || x == -z || y == z || y == -z || x + y == 0) return true;
  // the four L4 = 0 conics (drafts/ag_research_notes.md N5)
  long c1 = y*y - x*y - x*z - y*z - 10*z*z;
  long c2 = x*x - x*y - x*z - y*z - 10*z*z;
  long c3 = x*x - x*y + x*z + y*z - 10*z*z;
  long c4 = y*y - x*y + x*z + y*z - 10*z*z;
  return c1 == 0 || c2 == 0 || c3 == 0 || c4 == 0;
}

int main(int argc, char** argv) {
  long H = argc > 1 ? atol(argv[1]) : 1000;
  long z0 = argc > 2 ? atol(argv[2]) : 1, z1 = argc > 3 ? atol(argv[3]) : H;
  for (int i = 0; i < 64; ++i) sq64[(i * i) % 64] = true;
  for (int i = 0; i < 63; ++i) sq63[(i * i) % 63] = true;
  for (int i = 0; i < 65; ++i) sq65[(i * i) % 65] = true;
  for (int i = 0; i < 11; ++i) sq11[(i * i) % 11] = true;
  long long total = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+:total)
  for (long z = z0; z <= z1; ++z) {
    const i128 Z = z, Z2 = Z * Z, Z3 = Z2 * Z, Z4 = Z2 * Z2, Z5 = Z4 * Z;
    std::vector<long> hits;
    for (long x = -H; x <= H; ++x) {
      const i128 X = x, X2 = X * X, X3 = X2 * X, X4 = X2 * X2;
      for (long y = -H; y <= x; ++y) {
        const i128 Y = y, Y2 = Y * Y, Y3 = Y2 * Y, Y4 = Y2 * Y2;
        i128 E = X4 * Y2 - 2 * X3 * Y3 + X2 * Y4
               + Z2 * (10 * X3 * Y - 20 * X2 * Y2 + 10 * X * Y3)
               + Z4 * (45 * X2 - 10 * X * Y + 45 * Y2);
        i128 O = -4 * Z3 * (X3 + X2 * Y + X * Y2 + Y3) + 40 * Z5 * (X + Y);
        i128 q1 = E + O, q2 = E - O;
        if (!maybe_square(q1) || !maybe_square(q2)) continue;
        if (!is_square(q1) || !is_square(q2)) continue;
        if (degenerate(x, y, z)) continue;
        if (std::gcd(std::gcd(std::labs(x), std::labs(y)), z) != 1) continue;
        hits.push_back(x); hits.push_back(y);
      }
    }
#pragma omp critical
    {
      for (size_t i = 0; i < hits.size(); i += 2) printf("%ld %ld %ld\n", hits[i], hits[i + 1], z);
      fflush(stdout);
    }
    total += hits.size() / 2;
  }
  fprintf(stderr, "H = %ld, z in [%ld, %ld]: %lld nondegenerate-line hits\n", H, z0, z1, total);
  return 0;
}
