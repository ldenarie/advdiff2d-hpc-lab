// Mini-STREAM: measures the sustainable main-memory bandwidth of the machine,
// to give the benchmark results a hardware ceiling. See docs/benchmarking.md.
//
// Usage: stream [n=N] [reps=R]
//   n      number of doubles per array (default 2^25 = 256 MiB per array)
//   reps   timed repetitions per kernel; the best one is reported (default 20)
//
// Kernels (the same as in STREAM, McCalpin):
//   scale   a[i] = s * b[i]          1 read, 1 write
//   triad   a[i] = b[i] + s * c[i]   2 reads, 1 write
//
// `scale` has the same minimal memory traffic as one stencil update (read the
// old value, write the new one), so its element rate (`melem_s`, million
// elements per second) is the highest update rate, in MLUP/s, that a
// memory-bound stencil can reach on large grids: it is the ceiling.
// A plain copy is not used: compilers turn it into memcpy, which may use
// non-temporal stores and skip the write-allocate traffic a stencil pays.
//
// Bandwidth is reported twice:
//   stream_gbs    STREAM convention, counting only the bytes the code names
//                 (scale: 16 B per element, triad: 24 B);
//   traffic_gbs   including the write-allocate read of the destination
//                 (scale: 24 B, triad: 32 B), the convention used for the
//                 stencil's effective bandwidth in bench/plot.py.
//
// Output on stdout, CSV:
//   kernel,n,threads,best_s,median_s,stream_gbs,traffic_gbs,melem_s

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <format>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

long parse_arg(int argc, char** argv, std::string_view key, long default_value) {
    for (int k = 1; k < argc; ++k) {
        const std::string_view arg = argv[k];
        if (arg.starts_with(key) && arg.size() > key.size() && arg[key.size()] == '=') {
            return std::stol(std::string(arg.substr(key.size() + 1)));
        }
    }
    return default_value;
}

int thread_count() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

template <class Kernel>
std::vector<double> time_kernel(Kernel&& kernel, int reps) {
    kernel();  // warm-up: page faults, frequency ramp-up
    std::vector<double> times;
    times.reserve(static_cast<std::size_t>(reps));
    for (int r = 0; r < reps; ++r) {
        const auto t0 = std::chrono::steady_clock::now();
        kernel();
        const auto t1 = std::chrono::steady_clock::now();
        times.push_back(std::chrono::duration<double>(t1 - t0).count());
    }
    std::sort(times.begin(), times.end());
    return times;
}

void report(std::string_view name, long n, const std::vector<double>& sorted_times,
            double stream_bytes, double traffic_bytes) {
    const double best = sorted_times.front();
    const double median = sorted_times[sorted_times.size() / 2];
    const double elements = static_cast<double>(n);
    std::cout << std::format("{},{},{},{:.6f},{:.6f},{:.2f},{:.2f},{:.1f}\n", name, n,
                             thread_count(), best, median,
                             stream_bytes * elements / best / 1e9,
                             traffic_bytes * elements / best / 1e9,
                             elements / best / 1e6);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const long n = parse_arg(argc, argv, "n", 1L << 25);
        const int reps = static_cast<int>(parse_arg(argc, argv, "reps", 20));
        if (n <= 0 || reps <= 0) {
            throw std::invalid_argument("n and reps must be > 0");
        }

        // Not std::vector: its value-initialization would touch every page from
        // one thread. Each thread first-touches the part it will use.
        auto a = std::make_unique_for_overwrite<double[]>(static_cast<std::size_t>(n));
        auto b = std::make_unique_for_overwrite<double[]>(static_cast<std::size_t>(n));
        auto c = std::make_unique_for_overwrite<double[]>(static_cast<std::size_t>(n));
        double* const pa = a.get();
        double* const pb = b.get();
        double* const pc = c.get();

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (long i = 0; i < n; ++i) {
            pa[i] = 0.0;
            pb[i] = 1.0;
            pc[i] = 2.0;
        }

        const double s = 3.0;
        auto scale = [=] {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (long i = 0; i < n; ++i) {
                pa[i] = s * pb[i];
            }
        };
        auto triad = [=] {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (long i = 0; i < n; ++i) {
                pa[i] = pb[i] + s * pc[i];
            }
        };

        std::cout << "kernel,n,threads,best_s,median_s,stream_gbs,traffic_gbs,melem_s\n";
        report("scale", n, time_kernel(scale, reps), 16.0, 24.0);
        report("triad", n, time_kernel(triad, reps), 24.0, 32.0);

        // Use the result so that the compiler cannot discard the loops.
        if (pa[n / 2] != pb[n / 2] + s * pc[n / 2]) {
            throw std::runtime_error("unexpected result");
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& e) {
        std::cerr << std::format("error: {}\n", e.what());
        return EXIT_FAILURE;
    }
}
