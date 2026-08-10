// GPU measurement discipline.
//
// ---------------------------------------------------------------------------
// The thing that invalidates naive GPU benchmarks
// ---------------------------------------------------------------------------
//
// This RTX 4060 Laptop idles at **210 MHz against a 3105 MHz maximum** with
// persistence mode disabled. A kernel timed on a cold device therefore runs
// roughly 15x slow, and -- worse than being wrong -- it is wrong *silently* and
// *reproducibly*, so repeating the measurement confirms it.
//
// So every measurement here:
//
//   1. warms the device until the SM clock stops rising (NVML, polled between
//      launches), and
//   2. prints the clock it actually achieved, next to the result.
//
// A GPU number in this project without a clock beside it should not be
// believed. That is why Report() takes the clock rather than making it
// optional.
//
// ---------------------------------------------------------------------------
// Rooflines come from the device, not from a table
// ---------------------------------------------------------------------------
//
// SM count, clock and memory bus are queried through cudaGetDeviceProperties so
// the ceilings follow the hardware rather than a constant someone typed. The
// one number that must be supplied by hand is INT32 lanes per SM: Ada devotes
// 128 FP32 units per SM of which 64 also do INT32, and no runtime query reports
// that.
//
// ---------------------------------------------------------------------------
// Nsight Compute
// ---------------------------------------------------------------------------
//
// `ncu` is present but blocked: ERR_NVGPUCTRPERM, because driver 550 defaults
// to NVreg_RestrictProfilingToAdminUsers=1. Either run it under sudo, or to
// enable it for everyone:
//
//   echo 'options nvidia NVreg_RestrictProfilingToAdminUsers=0' |
//     sudo tee /etc/modprobe.d/nvidia-profiling.conf
//   sudo update-initramfs -u && reboot
//
// Nothing here needs it -- CUDA events and NVML are unprivileged -- but warp
// stall reasons and memory-pipe utilisation are only reachable that way.
#pragma once

#include <cuda_runtime.h>
#include <nvml.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace pte::gpubench {

#define PTE_CUDA_OK(call)                                                      \
  do {                                                                         \
    const cudaError_t err_ = (call);                                           \
    if (err_ != cudaSuccess) {                                                 \
      std::fprintf(stderr, "%s:%d: %s -> %s\n", __FILE__, __LINE__, #call,     \
                   cudaGetErrorString(err_));                                  \
      std::exit(1);                                                            \
    }                                                                          \
  } while (0)

// ---------------------------------------------------------------------------
// Device facts and rooflines
// ---------------------------------------------------------------------------

struct Device {
  cudaDeviceProp prop{};
  bool nvml = false;

  // Ada (sm_89): 128 FP32 lanes per SM, of which 64 also execute INT32. This is
  // the *nominal* figure and P1 shows it is a lower bound, not a ceiling:
  // mul.lo.u32 measures at 99.4% of it (so IMAD really is 64/SM/clk), but
  // add.u32 measures at 177%, meaning simple integer adds issue on more than the
  // INT32-designated lanes. A "%nominal" above 100 is therefore information
  // about the hardware, not a broken benchmark -- unlike the 1231% a deleted
  // loop produced, which is why P1 emits its operations as asm volatile.
  static constexpr double kInt32LanesPerSM = 64.0;

  static Device Query() {
    Device d;
    PTE_CUDA_OK(cudaGetDeviceProperties(&d.prop, 0));
    d.nvml = nvmlInit_v2() == NVML_SUCCESS;
    return d;
  }

  double MaxClockGHz() const { return prop.clockRate / 1e6; }

  // Peak 32-bit integer ops/s at a given achieved clock.
  double Int32Roofline(double ghz) const {
    return prop.multiProcessorCount * kInt32LanesPerSM * ghz * 1e9;
  }

  // Peak device memory bandwidth, bytes/s. memoryClockRate is in kHz and is the
  // *data* rate already; GDDR6 transfers on both edges, hence the 2.
  double MemoryRoofline() const {
    return 2.0 * prop.memoryClockRate * 1e3 * (prop.memoryBusWidth / 8.0);
  }

  // Current SM clock in GHz, or 0 if NVML is unavailable.
  double CurrentGHz() const {
    if (!nvml)
      return 0;
    nvmlDevice_t dev;
    if (nvmlDeviceGetHandleByIndex_v2(0, &dev) != NVML_SUCCESS)
      return 0;
    unsigned int mhz = 0;
    if (nvmlDeviceGetClockInfo(dev, NVML_CLOCK_SM, &mhz) != NVML_SUCCESS)
      return 0;
    return mhz / 1000.0;
  }

  void Print() const {
    std::printf("device: %s, sm_%d%d, %d SMs, %.0f MB, "
                "max SM clock %.0f MHz, bus %d-bit\n",
                prop.name, prop.major, prop.minor, prop.multiProcessorCount,
                prop.totalGlobalMem / 1048576.0, prop.clockRate / 1e3,
                prop.memoryBusWidth);
    std::printf("rooflines: INT32 %.1f Gop/s at max clock, memory %.0f GB/s\n",
                Int32Roofline(MaxClockGHz()) / 1e9, MemoryRoofline() / 1e9);
    if (!nvml)
      std::printf("WARNING: NVML unavailable -- clocks cannot be reported, so "
                  "every number below is suspect\n");
  }
};

// ---------------------------------------------------------------------------
// Timing
// ---------------------------------------------------------------------------

struct Measurement {
  double seconds = 0; // best of the measured iterations
  double ghz = 0;     // SM clock achieved during the measurement
  double ghz_start = 0;
};

// Runs `launch` until the clock plateaus, then measures it.
//
// The warm-up condition is the clock, not an iteration count: on a cold device
// the ramp takes tens of milliseconds, and on a hot one it never happens. Both
// are handled by watching the clock and stopping when it stops rising.
template <typename F>
Measurement Time(const Device &dev, int iters, F &&launch) {
  cudaEvent_t a, b;
  PTE_CUDA_OK(cudaEventCreate(&a));
  PTE_CUDA_OK(cudaEventCreate(&b));

  const double ghz_cold = dev.CurrentGHz();

  // Warm up: launch until the SM clock stops climbing, capped so a device that
  // reports nothing still terminates.
  double last = 0;
  for (int i = 0; i < 200; ++i) {
    launch();
    PTE_CUDA_OK(cudaDeviceSynchronize());
    const double now = dev.CurrentGHz();
    if (now > 0 && now <= last * 1.01 && i >= 3)
      break;
    last = now;
  }

  double best = 1e300;
  double ghz_sum = 0;
  int ghz_n = 0;
  for (int i = 0; i < iters; ++i) {
    PTE_CUDA_OK(cudaEventRecord(a));
    launch();
    PTE_CUDA_OK(cudaEventRecord(b));
    PTE_CUDA_OK(cudaEventSynchronize(b));
    float ms = 0;
    PTE_CUDA_OK(cudaEventElapsedTime(&ms, a, b));
    best = std::min(best, ms / 1e3);
    const double g = dev.CurrentGHz();
    if (g > 0) {
      ghz_sum += g;
      ++ghz_n;
    }
  }
  PTE_CUDA_OK(cudaEventDestroy(a));
  PTE_CUDA_OK(cudaEventDestroy(b));
  return {best, ghz_n ? ghz_sum / ghz_n : 0.0, ghz_cold};
}

// ---------------------------------------------------------------------------
// Reporting
// ---------------------------------------------------------------------------

class Table {
public:
  Table(std::string title, std::string unit)
      : title_(std::move(title)), unit_(std::move(unit)) {
    std::printf("\n=== %s ===\n", title_.c_str());
    std::printf("  %-38s %12s %12s %8s %9s\n", "case", "time", unit_.c_str(),
                "GHz", "%nominal");
  }

  void Row(const std::string &name, const Measurement &m, double units,
           double roofline_rate) const {
    const double rate = units / m.seconds;
    std::printf("  %-38s %9.3f ms %12.4g %8.2f", name.c_str(),
                m.seconds * 1e3, rate, m.ghz);
    if (roofline_rate > 0) {
      const double pct = 100.0 * rate / roofline_rate;
      std::printf(" %8.1f%%%s", pct, pct > 105.0 ? "  > nominal" : "");
    }
    std::printf("\n");
  }

  void Note(const std::string &s) const {
    std::printf("  %s\n", s.c_str());
  }

private:
  std::string title_;
  std::string unit_;
};

} // namespace pte::gpubench
