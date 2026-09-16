#pragma once

#include <stdexcept>
#include <string>
#include <functional>
#include <cstddef>

namespace flatland {

// The failure that ended a parallel_for, carrying which item produced it so the
// caller can phrase the message in its own terms ("timestep 7", "view 7", ...).
struct ParallelError : std::runtime_error {
    size_t index;
    ParallelError(size_t i, const std::string& what) : std::runtime_error(what), index(i) {}
};

// Run body(k, worker) for every k in [0, n), across at most n_workers threads.
//
// `worker` is a stable index in [0, n_workers), so a caller can keep per-worker
// scratch — a Renderer, a field buffer — in a plain vector indexed by it,
// without thread_local storage and without locking.
//
// If a body throws, the remaining work is abandoned and the LOWEST-indexed
// failure is rethrown as a ParallelError once every worker has joined. Reporting
// whichever thread happened to fail first would make the same input produce
// different messages from run to run.
//
// Work is handed out by an atomic counter rather than partitioned up front, so
// views of wildly different cost still spread evenly across the workers.
void parallel_for(size_t n, unsigned n_workers,
                  const std::function<void(size_t k, unsigned worker)>& body);

// How many workers to use for n items given a requested thread count
// (0 meaning "one per core"). Never exceeds n: spawning threads with nothing
// to do only costs.
unsigned choose_workers(size_t n, unsigned requested);

} // namespace flatland
