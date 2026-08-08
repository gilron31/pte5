// Host side of the E_4 meet: validates the batch, stages it through pinned
// memory, moves it to the device, launches, and copies the result back.
#include "gpu/meet/meet.h"

#include "gpu/core/device_append_buffer.cuh"
#include "gpu/core/device_buffer.cuh"
#include "gpu/core/host_buffer.cuh"
#include "gpu/meet/meet_impl.cuh"

#include <algorithm>
#include <stdexcept>

namespace pte::gpu {
namespace {

// Why is this a .cu ? it should be a .cc

// Staging a pageable input through pinned memory before the H2D copy lets the
// driver DMA it directly instead of bouncing it through its own pinned buffer.
// See host_buffer.cuh, and apps/hybrid_demo.cc for the measured difference.
template <typename T>
void CopyViaPinned(const T *src, std::size_t count, DeviceBuffer<T> &dst) {
  if (count == 0)
    return;
  HostBuffer<T> staging(count);
  std::copy(src, src + count, staging.begin());
  dst.CopyFromHost(staging.data(), count);
}

} // namespace

index_t meet_batch(const nu_t *nu, index_t nu_count, const index_t *offsets,
                   index_t num_norms, MeetHit *hits, index_t max_hits) {
  if (num_norms < 0 || max_hits < 0 || nu_count < 0)
    throw std::invalid_argument("meet_batch: negative size");
  if (num_norms == 0)
    return 0;
  if (nu == nullptr || offsets == nullptr)
    throw std::invalid_argument("meet_batch: null input");
  if (max_hits > 0 && hits == nullptr)
    throw std::invalid_argument("meet_batch: null hits buffer");

  // The kernel indexes nu through offsets without any further bound, so a
  // malformed offsets array would be an out-of-range device read. Reject it
  // here instead.
  if (offsets[0] < 0 || offsets[num_norms] > nu_count)
    throw std::invalid_argument("meet_batch: offsets out of range");
  for (index_t i = 0; i < num_norms; ++i) {
    if (offsets[i] > offsets[i + 1])
      throw std::invalid_argument("meet_batch: offsets not non-decreasing");
  }

  const std::size_t nu_used = static_cast<std::size_t>(offsets[num_norms]);

  DeviceBuffer<nu_t> d_nu(nu_used);
  DeviceBuffer<index_t> d_offsets(static_cast<std::size_t>(num_norms) + 1);
  DeviceAppendBuffer<MeetHit> d_hits(max_hits);

  CopyViaPinned(nu, nu_used, d_nu);
  CopyViaPinned(offsets, static_cast<std::size_t>(num_norms) + 1, d_offsets);

  LaunchE4MeetDirect(d_nu.get(), d_offsets.get(), num_norms, d_hits.view());

  // Reading the counter synchronizes, which is where an execution error would
  // surface. Only the prefix that was actually stored is copied back.
  const count_t found = d_hits.total();
  d_hits.CopyToHost(hits, d_hits.stored(found));
  return static_cast<index_t>(found);
}

} // namespace pte::gpu
