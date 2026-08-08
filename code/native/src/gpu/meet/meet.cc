// Host side of the E_4 meet: validates the batch, stages it through pinned
// memory, moves it to the device, launches, and copies the result back.
#include "gpu/meet/meet.h"

#include "gpu/core/device_append_buffer.h"
#include "gpu/core/device_buffer.h"
#include "gpu/core/host_buffer.h"
#include "gpu/meet/meet_impl.h"

#include <algorithm>
#include <stdexcept>

namespace pte::gpu {
namespace {

// Staging a pageable input through pinned memory before the H2D copy lets the
// driver DMA it directly instead of bouncing it through its own pinned buffer.
// See host_buffer.h, and apps/hybrid_demo.cc for the measured difference.
template <typename T>
void CopyViaPinned(const T *src, std::size_t count, DeviceBuffer<T> &dst) {
  if (count == 0)
    return;
  HostBuffer<T> staging(count);
  std::copy(src, src + count, staging.begin());
  dst.CopyFromHost(staging.data(), count);
}

} // namespace

index_t meet_batch(std::span<const nu_t> nu, std::span<const index_t> offsets,
                   std::span<MeetHit> hits) {
  if (offsets.size() < 2)
    return 0; // fewer than two offsets describes no norms

  // The kernel indexes nu through offsets without any further bound, so a
  // malformed offsets array would be an out-of-range device read. Reject it
  // here instead. Sizes themselves need no checking any more: a span cannot
  // misreport the extent of its own data.
  if (offsets.front() < 0 ||
      static_cast<std::size_t>(offsets.back()) > nu.size())
    throw std::invalid_argument("meet_batch: offsets out of range");
  for (std::size_t i = 0; i + 1 < offsets.size(); ++i) {
    if (offsets[i] > offsets[i + 1])
      throw std::invalid_argument("meet_batch: offsets not non-decreasing");
  }

  const auto nu_used = static_cast<std::size_t>(offsets.back());
  const auto max_hits = static_cast<index_t>(hits.size());

  DeviceBuffer<nu_t> d_nu(nu_used);
  DeviceBuffer<index_t> d_offsets(offsets.size());
  DeviceAppendBuffer<MeetHit> d_hits(max_hits);

  CopyViaPinned(nu.data(), nu_used, d_nu);
  CopyViaPinned(offsets.data(), offsets.size(), d_offsets);

  LaunchE4MeetDirect(d_nu, d_offsets, d_hits.view());

  // Reading the counter synchronizes, which is where an execution error would
  // surface. Only the prefix that was actually stored is copied back.
  const count_t found = d_hits.total();
  d_hits.CopyToHost(hits.data(), d_hits.stored(found));
  return static_cast<index_t>(found);
}

} // namespace pte::gpu
