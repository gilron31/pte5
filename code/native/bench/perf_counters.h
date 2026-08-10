// Hardware performance counters for a code region.
//
// Wraps perf_event_open(2) directly rather than shelling out to `perf stat`,
// and the reason is the whole point of this file: `perf stat` measures a
// *process*, but every interesting question here is about a *phase*. The walk
// is already decomposed by ablation (bench_walk_phases), so scoping counters to
// the same regions makes per-phase IPC fall out of a measurement we are taking
// anyway. A process-level tool cannot do that without a second run per phase.
//
// ---------------------------------------------------------------------------
// Availability
// ---------------------------------------------------------------------------
//
// Needs kernel.perf_event_paranoid <= 1 (Ubuntu ships 4, meaning "disallow all
// unprivileged"). If the syscall fails this class stays quiet and reports
// `available() == false` with a reason -- there is no second code path and no
// #ifdef, because a benchmark that only compiles on a tuned machine is a
// benchmark nobody runs.
//
//   sudo sysctl -w kernel.perf_event_paranoid=1
//   echo 'kernel.perf_event_paranoid=1' | sudo tee /etc/sysctl.d/60-perf.conf
//
// ---------------------------------------------------------------------------
// Why a counter *group*
// ---------------------------------------------------------------------------
//
// The five events are opened as one group with cycles as leader, so the kernel
// schedules them onto the PMU together and they all cover exactly the same
// instructions. Opened separately they could be multiplexed onto different time
// slices, and then `instructions / cycles` would be a ratio of two numbers
// measured over different windows -- an IPC that is quietly meaningless. Zen 4
// has six general-purpose counters, so five events need no multiplexing at all.
//
// Counting is user-space only (exclude_kernel), which keeps a page fault or a
// scheduler tick out of the numbers.
#pragma once

#include <asm/unistd.h>
#include <linux/perf_event.h>
#include <sys/ioctl.h>
#include <unistd.h>

#include <cerrno>
#include <cstdint>
#include <cstring>
#include <string>

namespace pte::bench {

class PerfCounters {
public:
  struct Reading {
    bool available = false;
    std::uint64_t cycles = 0;
    std::uint64_t instructions = 0;
    std::uint64_t branch_misses = 0;
    std::uint64_t cache_misses = 0;
    std::uint64_t l1d_load_misses = 0;

    double ipc() const {
      return cycles ? static_cast<double>(instructions) /
                          static_cast<double>(cycles)
                    : 0.0;
    }
    // Cycles per unit of whatever the caller is counting (norms, partitions).
    double per_unit(double units) const {
      return units > 0 ? static_cast<double>(cycles) / units : 0.0;
    }
  };

  PerfCounters() { Open(); }
  ~PerfCounters() { Close(); }
  PerfCounters(const PerfCounters &) = delete;
  PerfCounters &operator=(const PerfCounters &) = delete;

  bool available() const { return fds_[0] >= 0; }
  const std::string &reason() const { return reason_; }

  void Start() {
    if (!available())
      return;
    ioctl(fds_[0], PERF_EVENT_IOC_RESET, PERF_IOC_FLAG_GROUP);
    ioctl(fds_[0], PERF_EVENT_IOC_ENABLE, PERF_IOC_FLAG_GROUP);
  }

  Reading Stop() {
    Reading r;
    if (!available())
      return r;
    ioctl(fds_[0], PERF_EVENT_IOC_DISABLE, PERF_IOC_FLAG_GROUP);

    // PERF_FORMAT_GROUP layout: { u64 nr; u64 values[nr]; } in the order the
    // events were added to the group.
    struct {
      std::uint64_t nr;
      std::uint64_t values[kEvents];
    } buf{};
    const ssize_t n = read(fds_[0], &buf, sizeof(buf));
    if (n < 0 || buf.nr != kEvents)
      return r;

    r.available = true;
    r.cycles = buf.values[0];
    r.instructions = buf.values[1];
    r.branch_misses = buf.values[2];
    r.cache_misses = buf.values[3];
    r.l1d_load_misses = buf.values[4];
    return r;
  }

private:
  static constexpr int kEvents = 5;

  static long OpenOne(perf_event_attr &attr, int group_fd) {
    return syscall(__NR_perf_event_open, &attr, /*pid=*/0, /*cpu=*/-1, group_fd,
                   /*flags=*/0UL);
  }

  static perf_event_attr Attr(std::uint32_t type, std::uint64_t config,
                              bool leader) {
    perf_event_attr a{};
    a.size = sizeof(a);
    a.type = type;
    a.config = config;
    a.disabled = leader ? 1 : 0;
    a.exclude_kernel = 1;
    a.exclude_hv = 1;
    a.inherit = 0; // this thread only; each bench thread opens its own
    if (leader)
      a.read_format = PERF_FORMAT_GROUP;
    return a;
  }

  static constexpr std::uint64_t kL1dReadMiss =
      PERF_COUNT_HW_CACHE_L1D |
      (static_cast<std::uint64_t>(PERF_COUNT_HW_CACHE_OP_READ) << 8) |
      (static_cast<std::uint64_t>(PERF_COUNT_HW_CACHE_RESULT_MISS) << 16);

  void Open() {
    struct Spec {
      std::uint32_t type;
      std::uint64_t config;
    };
    const Spec spec[kEvents] = {
        {PERF_TYPE_HARDWARE, PERF_COUNT_HW_CPU_CYCLES},
        {PERF_TYPE_HARDWARE, PERF_COUNT_HW_INSTRUCTIONS},
        {PERF_TYPE_HARDWARE, PERF_COUNT_HW_BRANCH_MISSES},
        {PERF_TYPE_HARDWARE, PERF_COUNT_HW_CACHE_MISSES},
        {PERF_TYPE_HW_CACHE, kL1dReadMiss},
    };

    for (int i = 0; i < kEvents; ++i) {
      perf_event_attr a = Attr(spec[i].type, spec[i].config, i == 0);
      const long fd = OpenOne(a, i == 0 ? -1 : fds_[0]);
      if (fd < 0) {
        reason_ = std::string("perf_event_open failed on event ") +
                  std::to_string(i) + ": " + std::strerror(errno) +
                  (errno == EACCES || errno == EPERM
                       ? " (try: sudo sysctl -w kernel.perf_event_paranoid=1)"
                       : "");
        Close();
        return;
      }
      fds_[i] = static_cast<int>(fd);
    }
  }

  void Close() {
    for (int &fd : fds_) {
      if (fd >= 0)
        close(fd);
      fd = -1;
    }
  }

  int fds_[kEvents] = {-1, -1, -1, -1, -1};
  std::string reason_ = "not opened";
};

} // namespace pte::bench
