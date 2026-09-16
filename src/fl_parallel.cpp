// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Interfluo
//
// FlatLand is dual-licensed: GNU AGPL v3 (see LICENSE) or a commercial licence
// for closed-source or hosted use (see COMMERCIAL-LICENSE.md).

#include "fl_parallel.hpp"

#include <algorithm>
#include <atomic>
#include <limits>
#include <mutex>
#include <system_error>
#include <thread>
#include <vector>

namespace flatland {

unsigned choose_workers(size_t n, unsigned requested) {
    if (n == 0) return 1;
    if (requested == 0) {
        const unsigned hw = std::thread::hardware_concurrency();
        requested = hw ? hw : 1;
    }
    const size_t capped = std::min<size_t>(requested, n);
    return (unsigned)std::max<size_t>(1, capped);
}

void parallel_for(size_t n, unsigned n_workers,
                  const std::function<void(size_t, unsigned)>& body) {
    if (n == 0) return;
    n_workers = std::max(1u, (unsigned)std::min<size_t>(n_workers, n));

    std::atomic<size_t> next{0};
    std::mutex err_mtx;
    std::string err_msg;
    size_t err_k = std::numeric_limits<size_t>::max();

    auto run = [&](unsigned w) {
        size_t k;
        while ((k = next.fetch_add(1)) < n) {
            try {
                body(k, w);
            } catch (const std::exception& e) {
                std::lock_guard<std::mutex> lk(err_mtx);
                if (k < err_k) { err_k = k; err_msg = e.what(); }
                next.store(n);              // tell the other workers to stop
                return;
            } catch (...) {
                std::lock_guard<std::mutex> lk(err_mtx);
                if (k < err_k) { err_k = k; err_msg = "unknown error"; }
                next.store(n);
                return;
            }
        }
    };

    if (n_workers <= 1) {
        run(0);
    } else {
        std::vector<std::thread> pool;
        pool.reserve(n_workers);
        try {
            for (unsigned w = 0; w < n_workers; ++w) pool.emplace_back(run, w);
        } catch (const std::system_error&) {
            // The OS refused a thread. Stop the work and join what did start:
            // letting this unwind through a vector still holding joinable
            // threads calls std::terminate.
            next.store(n);
            for (auto& th : pool) if (th.joinable()) th.join();
            throw std::runtime_error("could not start " + std::to_string(n_workers) +
                                     " worker threads; request fewer");
        }
        for (auto& th : pool) th.join();
    }

    if (err_k != std::numeric_limits<size_t>::max())
        throw ParallelError(err_k, err_msg);
}

} // namespace flatland
