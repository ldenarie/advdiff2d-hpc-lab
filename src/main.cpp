// Driver: sets up a Gaussian tracer, runs the scheme, benchmark the performances
//         over several runs (data is only set once)
//         u = 1.0 (velocity in the X axis)
//         v = 0.5 (velocity in the Y axis)
//         kappa = 1e-4 (diffusion coefficient)
//         s = 0.9 (stability number)
//
//         The field is allocated and initialized once; repetitions continue the same simulation.
//         Warm-up repetitions absorb first-touch page faults and cache effects of the first steps.
//
// Usage: advdiff2d [key=value ...]
//   kernel     used kernel (weighted|naive|incremental) (default incremental)
//   nx, ny     grid size                                (default 512 x 512)
//   steps      number of time steps in a run            (default 200)
//   reps       number of run to benchmark               (default to 3)
//   warmup     number of warmup run                     (default to 1)
//   output     output format (text|csv)                 (default to text)

#include <advdiff/field.hpp>
#include <advdiff/kernel.hpp>
#include <advdiff/problem.hpp>

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <format>
#include <functional>
#include <iostream>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

namespace {

using Options = std::map<std::string, std::string, std::less<>>;

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

template <class T>
T get(const Options& o, std::string_view key, std::optional<T> default_value = std::nullopt) {
    const auto it = o.find(key);
    if (it == o.end()) {
        if (default_value.has_value()) {
            return default_value.value();
        }
        throw std::invalid_argument(std::format("Missing argument: {}", key));
    }
    if constexpr (std::is_same_v<T, int>) {
        return std::stoi(it->second);
    } else if constexpr (std::is_same_v<T, std::string>) {
        return it->second;
    } else {
        static_assert(!sizeof(T), "get<T>: unsupported option type");
    }
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

struct BenchmarkConfiguration {
    std::string kernel;
    int steps;
    int warmup;
    int reps;
};

struct RunParameters {
    advdiff::Physics p;
    advdiff::Grid g;
    double dt;
    double s;
    int nb_steps;
};

struct FieldData {
    double mass_rel_err;
    double min_init;
    double min_final;
    double max_init;
    double max_final;
};

struct RunResults {
    int run_number;
    double run_duration;
    double mlups;
    FieldData field_data;
};

template <class OUTPUT_RECORDER>
void run_configuration(int nx, int ny, const BenchmarkConfiguration& conf, OUTPUT_RECORDER&& record_output) {
    RunParameters params{
        .p = {
            .u = 1.,
            .v = .5,
            .kappa = 1e-4
        },
        .g = {
            .nx = nx,
            .ny = ny
        },
        .dt = 1.,
        .s = .9,
        .nb_steps = conf.steps
    };

    // dt from the target stability number
    params.dt = params.s / advdiff::stability_number(params.g, params.p, 1.);

    advdiff::Field c(params.g.nx, params.g.ny);
    advdiff::Field c_new(params.g.nx, params.g.ny);
    init_gaussian(c, params.g, 0.5 * params.g.lx, 0.5 * params.g.ly, 0.05 * params.g.lx);

    auto advance = [&c, &c_new, &params, &conf](int nb_steps) {
        auto kernel_function = [](std::string_view kernel) {
            if (kernel =="naive") {
                return advdiff::step_naive;
            } else if (kernel == "weighted") {
                return advdiff::step_weighted;
            }
            return advdiff::step_incremental;
        } (conf.kernel);
        for (int n = 0; n < nb_steps; ++n) {
            advdiff::apply_periodic_bc(c);
            kernel_function(c, c_new, params.g, params.p, params.dt);
            std::swap(c, c_new); // O(1): swaps the underlying buffers
        }
    };

    // Run warmup rounds
    for (int w = 0; w < conf.warmup; ++w) {
        advance(params.nb_steps);
    }

    // Record each repetition

    for (int r = 0; r < conf.reps; ++r) {
        const auto mass0 = advdiff::total_mass(c, params.g);
        const auto min0 = advdiff::field_min(c);
        const auto max0 = advdiff::field_max(c);

        const auto t0 = std::chrono::steady_clock::now();

        advance(params.nb_steps);

        const auto t1 = std::chrono::steady_clock::now();
        const auto run_duration = std::chrono::duration<double>(t1 - t0).count();

        const auto mass = advdiff::total_mass(c, params.g);
        const auto final_min_value = advdiff::field_min(c);
        const auto final_max_value = advdiff::field_max(c);
        const auto mlups = static_cast<double>(params.g.nx) * params.g.ny * params.nb_steps / run_duration / 1e6;

        const RunResults results{
            .run_number = r,
            .run_duration = run_duration,
            .mlups = mlups,
            .field_data = {
                .mass_rel_err = std::abs(mass - mass0) / mass0,
                .min_init = min0,
                .min_final = final_min_value,
                .max_init = max0,
                .max_final = final_max_value
            }
        };

        record_output(conf.kernel, params, results);
    }
}

void write_report(std::string_view kernel, const RunParameters& params, const RunResults& results, std::ostream& out) {
    const advdiff::Physics& p = params.p;
    const advdiff::Grid& g = params.g;
    out << "Simulation parameters:\n";
    out << std::format("u: {:.3f}, v: {:.3f}, kappa: {:.6e}\n", p.u, p.v, p.kappa);
    out << std::format("nx: {}, ny: {}\n", g.nx, g.ny);
    out << std::format("dt: {:.6e}, s: {:.3f}\n", params.dt, params.s);
    out << std::format("nb_steps: {}\n", params.nb_steps);
    out << std::format("kernel: {}\n", kernel);
    out << std::format("rep: {}\n", results.run_number);
    out << std::format("time_s: {:.6f}\n", results.run_duration);
    out << std::format("mlups: {:.2f}\n", results.mlups);
    out << std::format("mass_rel_err: {:.3e}\n", results.field_data.mass_rel_err);
    out << std::format("field_min: {:.6e} (init: {:.6e})\n", results.field_data.min_final, results.field_data.min_init);
    out << std::format("field_max: {:.6e} (init: {:.6e})\n", results.field_data.max_final, results.field_data.max_init);
    out << std::endl;
}

void record_output_csv(std::string_view kernel, const RunParameters& params, const RunResults& results) {
    write_report(kernel, params, results, std::cerr);

    std::cout << std::format("{},{},{},{},{},{:.6f},{:.2f},{:.3e}\n",
        kernel,
        params.g.nx,
        params.g.ny,
        params.nb_steps,
        results.run_number,
        results.run_duration,
        results.mlups,
        results.field_data.mass_rel_err);
    std::cout << std::flush;
}

void record_output_text(std::string_view kernel, const RunParameters& params, const RunResults& results) {
    write_report(kernel, params, results, std::cout);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options opts = parse_args(argc, argv);

        const std::string kernel = get<std::string>(opts, "kernel", "incremental");
        if (kernel != "weighted" and kernel != "naive" and kernel != "incremental") {
            throw std::invalid_argument(std::format("Unknown kernel. Got: {}, Expected: weighted|naive|incremental", kernel));
        }

        const int nx = get<int>(opts, "nx", 512);
        if (nx <= 0) {
            throw std::invalid_argument(std::format("nx must be > 0. Got: {}", nx));
        }

        const int ny = get<int>(opts, "ny", 512);
        if (ny <= 0) {
            throw std::invalid_argument(std::format("ny must be > 0. Got: {}", ny));
        }

        const int steps = get<int>(opts, "steps", 200);
        if (steps <= 0) {
            throw std::invalid_argument(std::format("steps must be > 0. Got: {}", steps));
        }

        const int reps = get<int>(opts, "reps", 3);
        if (reps <= 0) {
            throw std::invalid_argument(std::format("reps must be > 0. Got: {}", reps));
        }

        const int warmup = get<int>(opts, "warmup", 1);
        if (warmup < 0) {
            throw std::invalid_argument(std::format("warmup must be >= 0. Got: {}", warmup));
        }

        const std::string output = get<std::string>(opts, "output", "text");
        if (output != "text" && output != "csv") {
            throw std::invalid_argument(std::format("output argument invalid. Got: {}, Expected: text | csv", output));
        }

        if (output == "csv") {
            // write headers
            std::cout << "kernel,nx,ny,steps,rep,time_s,mlups,mass_rel_err" << std::endl;
        }
        auto record_output = output == "csv" ? record_output_csv : record_output_text;

        const BenchmarkConfiguration conf{
            .kernel = kernel,
            .steps = steps,
            .warmup = warmup,
            .reps = reps
        };
        run_configuration(nx, ny, conf, record_output);

        return EXIT_SUCCESS;
    } catch (const std::exception& e) {
        std::cerr << std::format("error: {}\n", e.what());
        return EXIT_FAILURE;
    }
}
