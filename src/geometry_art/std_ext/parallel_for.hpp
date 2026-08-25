#ifndef GEOMETRY_ART_STD_EXT_PARALLEL_FOR_HPP_
#define GEOMETRY_ART_STD_EXT_PARALLEL_FOR_HPP_

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <thread>
#include <vector>

namespace geometry_art::std_ext {

namespace detail {

template<std::invocable<size_t> Body>
void run_range(const Body& body, size_t begin, size_t end) {
    for (size_t index = begin; index < end; ++index) {
        body(index);
    }
}

inline size_t configured_worker_count() {
    const char* override = std::getenv("GEOMETRY_ART_THREADS");

    if (override == nullptr) {
        return std::thread::hardware_concurrency();
    }

    return static_cast<size_t>(std::atoi(override));
}

inline size_t worker_count(size_t count) {
    static const size_t available = configured_worker_count();
    return std::min(count, std::max<size_t>(1, available));
}

} // namespace detail

template<std::invocable<size_t> Body>
void parallel_for(size_t count, const Body& body) {
    size_t workers = detail::worker_count(count);

    if (workers <= 1) {
        detail::run_range(body, 0, count);
        return;
    }

    size_t chunk = (count + workers - 1) / workers;
    std::vector<std::exception_ptr> failures(workers);
    std::vector<std::thread> threads;
    threads.reserve(workers - 1);

    auto run_chunk = [&](size_t worker) {
        size_t begin = std::min(count, worker * chunk);

        try {
            detail::run_range(body, begin, std::min(count, begin + chunk));
        } catch (...) {
            failures[worker] = std::current_exception();
        }
    };

    for (size_t worker = 1; worker < workers; ++worker) {
        threads.emplace_back(run_chunk, worker);
    }

    run_chunk(0);

    for (std::thread& thread : threads) {
        thread.join();
    }

    for (const std::exception_ptr& failure : failures) {
        if (failure) {
            std::rethrow_exception(failure);
        }
    }
}

} // namespace geometry_art::std_ext

#endif //GEOMETRY_ART_STD_EXT_PARALLEL_FOR_HPP_
