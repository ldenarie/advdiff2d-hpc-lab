// Driver: sets up a Gaussian tracer, runs the scheme, prints diagnostics and timing.
//
// Usage: advdiff2d [key=value ...]
//   nx, ny     grid size                    (default 512 x 512)
//   steps      number of time steps         (default 200)
//   u, v       velocity                     (default 1.0, 0.5)
//   kappa      diffusion coefficient        (default 1e-4)
//   s          target stability number S    (default 0.9, must be <= 1)

#include "field.hpp"
#include "kernel.hpp"
#include "problem.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>

namespace {

using Options = std::map<std::string, std::string>;

Options parse_args(int argc, char** argv) {
    Options opts;
    for (int k = 1; k < argc; ++k) {
        const std::string arg = argv[k];
        const auto eq = arg.find('=');
        if (eq == std::string::npos) {
            throw std::invalid_argument("expected key=value, got '" + arg + "'");
        }
        opts[arg.substr(0, eq)] = arg.substr(eq + 1);
    }
    return opts;
}

double get(const Options& o, const std::string& key, double fallback) {
    const auto it = o.find(key);
    return it == o.end() ? fallback : std::stod(it->second);
}

void init_gaussian(advdiff::Field& c, const advdiff::Grid& g,
                   double x0, double y0, double sigma) {
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            const double rx = g.x(i) - x0;
            const double ry = g.y(j) - y0;
            c(i, j) = std::exp(-(rx * rx + ry * ry) / (2.0 * sigma * sigma));
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options opts = parse_args(argc, argv);

        advdiff::Grid g;
        g.nx = static_cast<int>(get(opts, "nx", 512));
        g.ny = static_cast<int>(get(opts, "ny", 512));
        const int steps = static_cast<int>(get(opts, "steps", 200));

        advdiff::Physics p;
        p.u = get(opts, "u", 1.0);
        p.v = get(opts, "v", 0.5);
        p.kappa = get(opts, "kappa", 1e-4);

        // dt from the target stability number: S is linear in dt.
        const double s_target = get(opts, "s", 0.9);
        const double dt = s_target / advdiff::stability_number(g, p, 1.0);

        advdiff::Field c(g.nx, g.ny);
        advdiff::Field c_new(g.nx, g.ny);
        init_gaussian(c, g, 0.5 * g.lx, 0.5 * g.ly, 0.05 * g.lx);

        const double mass0 = advdiff::total_mass(c, g);
        const double min0 = advdiff::field_min(c);
        const double max0 = advdiff::field_max(c);

        const auto t0 = std::chrono::steady_clock::now();
        for (int n = 0; n < steps; ++n) {
            advdiff::apply_periodic_bc(c);
            advdiff::step_seq(c, c_new, g, p, dt);
            std::swap(c, c_new);  // O(1): swaps the underlying buffers
        }
        const auto t1 = std::chrono::steady_clock::now();
        const double seconds = std::chrono::duration<double>(t1 - t0).count();

        const double mass = advdiff::total_mass(c, g);
        const double mlups = static_cast<double>(g.nx) * g.ny * steps / seconds / 1e6;

        std::printf("nx=%d ny=%d steps=%d dt=%.6e S=%.3f\n", g.nx, g.ny, steps, dt, s_target);
        std::printf("time_s=%.6f mlups=%.2f\n", seconds, mlups);
        std::printf("mass_rel_err=%.3e min=%.6e (init %.6e) max=%.6e (init %.6e)\n",
                    std::abs(mass - mass0) / mass0, advdiff::field_min(c), min0,
                    advdiff::field_max(c), max0);
        return EXIT_SUCCESS;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return EXIT_FAILURE;
    }
}
