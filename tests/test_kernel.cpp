// Verification tests. See docs/model.md, section 5, for the expected properties.
//
// TODO(block 2): implement each test and remove its SKIP().

#include <advdiff/field.hpp>
#include <advdiff/kernel.hpp>
#include <advdiff/problem.hpp>

#include <advdiff/details/kernel.hpp>

#include <catch2/catch_get_random_seed.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <random>

using advdiff::Field;
using advdiff::Grid;
using advdiff::Physics;

constexpr double ghost_cell_test_grid_value(int i, int j) {
    return 100 * i + j;
}
constexpr double ghost_cell_sentinel = -100;

TEST_CASE("Ghost cells copy the opposite interior edges", "[bc]") {
    const auto [nx, ny] = GENERATE(table<int, int>({
        {1, 1}, // Square: Limit case
        {2, 2}, // Square: Limit case
        {10, 10}, // Square: Normal case
        {2, 1}, // Rectangle: Limit case
        {1, 2}, // Rectangle: Limit case
        {10, 1}, // Rectangle: Limit case in y direction
        {1, 10}, // Rectangle: Limit case in x direction
        {10, 2}, // Rectangle: Limit case in y direction
        {2, 10}, // Rectangle: Limit case in x direction
        {4, 10}, // Rectangle: General case
        {10, 4} // Rectangle: General case
    }));
    CAPTURE(nx, ny);

    Field c(nx, ny);
    for (int j = 0; j < ny; ++j) {
        c(-1, j) = ghost_cell_sentinel;
        c(nx, j) = ghost_cell_sentinel;
        for (int i = 0; i < nx ; ++i) {
            c(i, j) = ghost_cell_test_grid_value(i, j);
        }
    }
    for (int i = 0; i < nx; ++i) {
        c(i, -1) = ghost_cell_sentinel;
        c(i, ny) = ghost_cell_sentinel;
    }

    advdiff::apply_periodic_bc(c);

    // Requires that the inner values did not change
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            REQUIRE(c(i, j) == ghost_cell_test_grid_value(i, j));
        }
    }

    // Requires ghost cell to have the value of the opposite cell in the grid
    for (int i = 0; i < nx; ++i) {
        REQUIRE(c(i, -1) == c(i, ny - 1));
        REQUIRE(c(i, ny) == c(i, 0));
    }
    for (int j = 0; j < ny; ++j) {
        REQUIRE(c(-1, j) == c(nx - 1, j));
        REQUIRE(c(nx, j) == c(0, j));
    }
}

TEST_CASE("Kernel function is correct", "[kernel]") {
    using advdiff::details::advdiff_kernel;

    advdiff::Physics p;
    p.u = GENERATE(-2., 0, 2.);
    p.v = GENERATE(-2., 0, 2.);
    p.kappa = GENERATE(0., 1e-3, 1., 2.);

    advdiff::Grid g;
    const auto dx = 1e-3;
    const auto dy = 1e-3;
    g.nx = 100;
    g.ny = 100;
    g.lx = g.nx * dx;
    g.ly = g.ny * dy;

    const auto dt = 1e-3;
    const auto cell = GENERATE(-2., 0, 2.);
    const auto west = GENERATE(-2., 0, 2.);
    const auto east = GENERATE(-2., 0, 2.);
    const auto south = GENERATE(-2., 0, 2.);
    const auto north = GENERATE(-2., 0, 2.);

    CAPTURE(p.u, p.v, p.kappa, cell, west, east, south, north);

    const double c_x = std::abs(p.u) * dt / dx;
    const double c_y = std::abs(p.v) * dt / dy;
    const double d_x = p.kappa * dt / (dx * dx);
    const double d_y = p.kappa * dt / (dy * dy);
    const double s = c_x + c_y + 2 * d_x + 2 * d_y;

    const double west_factor = p.u < 0 ? d_x : c_x + d_x;
    const double east_factor = p.u < 0 ? c_x + d_x : d_x;
    const double south_factor = p.v < 0 ? d_y : c_y + d_y;
    const double north_factor = p.v < 0 ? c_y + d_y : d_y;
    const double expected = (1 - s) * cell +
                            west_factor * west +
                            east_factor * east +
                            south_factor * south +
                            north_factor * north;

    REQUIRE_THAT(advdiff_kernel(cell, west, east, south, north, g, p, dt), Catch::Matchers::WithinAbs(expected, 1e-11));
}

TEST_CASE("step_seq for specific fields") {
    const auto rows = 3;
    const auto cols = 4;
    advdiff::Physics p;
    p.kappa = 1.;

    advdiff::Grid g;
    g.nx = cols;
    g.ny = rows;
    g.lx = cols; // dx = 1.
    g.ly = 2. * rows; // dy = 2.

    const auto dt = 1.;

    advdiff::Field c{g.nx, g.ny};
    c(0, 0) = 1.;
    c(1, 0) = -5.;
    c(2, 0) = 13.;
    c(3, 0) = 17.;
    c(0, 1) = -12.;
    c(1, 1) = 8.;
    c(2, 1) = -3.;
    c(3, 1) = -4.;
    c(0, 2) = -19.;
    c(1, 2) = 2.;
    c(2, 2) = 2.;
    c(3, 2) = -4.;

    advdiff::Field c_expected{g.nx, g.ny};

    SECTION("u > 0 and v > 0") {
        p.u = 1.;
        p.v = 3.;

        c_expected(0, 0) = -11.25;
        c_expected(1, 0) = 40.5;
        c_expected(2, 0) = -42.25;
        c_expected(3, 0) = -49.0;
        c_expected(0, 1) = 45.0;
        c_expected(1, 1) = -67.25;
        c_expected(2, 1) = 47.25;
        c_expected(3, 1) = 26.75;
        c_expected(0, 2) = 49.25;
        c_expected(1, 2) = -31.25;
        c_expected(2, 2) = -10.0;
        c_expected(3, 2) = -1.75;
    }

    SECTION("u > 0 and v = 0") {
        p.u = 1.;
        p.v = 0.;

        c_expected(0, 0) = 18.75;
        c_expected(1, 0) = 30.0;
        c_expected(2, 0) = -25.75;
        c_expected(3, 0) = -17.5;
        c_expected(0, 1) = 25.5;
        c_expected(1, 1) = -47.75;
        c_expected(2, 1) = 23.25;
        c_expected(3, 1) = -4.75;
        c_expected(0, 2) = 38.75;
        c_expected(1, 2) = -40.25;
        c_expected(2, 2) = -2.5;
        c_expected(3, 2) = -1.75;
    }

    SECTION("u > 0 and v < 0") {
        p.u = 1.;
        p.v = -3.;

        c_expected(0, 0) = -0.75;
        c_expected(1, 0) = 49.5;
        c_expected(2, 0) = -49.75;
        c_expected(3, 0) = -49.0;
        c_expected(0, 1) = 15.0;
        c_expected(1, 1) = -56.75;
        c_expected(2, 1) = 30.75;
        c_expected(3, 1) = -4.75;
        c_expected(0, 2) = 68.75;
        c_expected(1, 2) = -50.75;
        c_expected(2, 2) = 14.0;
        c_expected(3, 2) = 29.75;
    }

    SECTION("u = 0 and v > 0") {
        p.u = 0.;
        p.v = 3.;

        c_expected(0, 0) = -27.25;
        c_expected(1, 0) = 34.5;
        c_expected(2, 0) = -24.25;
        c_expected(3, 0) = -45.0;
        c_expected(0, 1) = 37.0;
        c_expected(1, 1) = -47.25;
        c_expected(2, 1) = 36.25;
        c_expected(3, 1) = 25.75;
        c_expected(0, 2) = 34.25;
        c_expected(1, 2) = -10.25;
        c_expected(2, 2) = -10.0;
        c_expected(3, 2) = -7.75;
    }

    SECTION("u = 0 and v = 0") {
        p.u = 0.;
        p.v = 0.;

        c_expected(0, 0) = 2.75;
        c_expected(1, 0) = 24.0;
        c_expected(2, 0) = -7.75;
        c_expected(3, 0) = -13.5;
        c_expected(0, 1) = 17.5;
        c_expected(1, 1) = -27.75;
        c_expected(2, 1) = 12.25;
        c_expected(3, 1) = -5.75;
        c_expected(0, 2) = 23.75;
        c_expected(1, 2) = -19.25;
        c_expected(2, 2) = -2.5;
        c_expected(3, 2) = -7.75;
    }

    SECTION("u = 0 and v < 0") {
        p.u = 0.;
        p.v = -3.;

        c_expected(0, 0) = -16.75;
        c_expected(1, 0) = 43.5;
        c_expected(2, 0) = -31.75;
        c_expected(3, 0) = -45.0;
        c_expected(0, 1) = 7.0;
        c_expected(1, 1) = -36.75;
        c_expected(2, 1) = 19.75;
        c_expected(3, 1) = -5.75;
        c_expected(0, 2) = 53.75;
        c_expected(1, 2) = -29.75;
        c_expected(2, 2) = 14.0;
        c_expected(3, 2) = 23.75;
    }

    SECTION("u < 0 and v > 0") {
        p.u = -1.;
        p.v = 3.;

        c_expected(0, 0) = -33.25;
        c_expected(1, 0) = 52.5;
        c_expected(2, 0) = -20.25;
        c_expected(3, 0) = -61.0;
        c_expected(0, 1) = 57.0;
        c_expected(1, 1) = -58.25;
        c_expected(2, 1) = 35.25;
        c_expected(3, 1) = 17.75;
        c_expected(0, 2) = 55.25;
        c_expected(1, 2) = -10.25;
        c_expected(2, 2) = -16.0;
        c_expected(3, 2) = -22.75;
    }

    SECTION("u < 0 and v = 0") {
        p.u = -1.;
        p.v = 0.;

        c_expected(0, 0) = -3.25;
        c_expected(1, 0) = 42.0;
        c_expected(2, 0) = -3.75;
        c_expected(3, 0) = -29.5;
        c_expected(0, 1) = 37.5;
        c_expected(1, 1) = -38.75;
        c_expected(2, 1) = 11.25;
        c_expected(3, 1) = -13.75;
        c_expected(0, 2) = 44.75;
        c_expected(1, 2) = -19.25;
        c_expected(2, 2) = -8.5;
        c_expected(3, 2) = -22.75;
    }

    SECTION("u < 0 and v < 0") {
        p.u = -1.;
        p.v = -3.;

        c_expected(0, 0) = -22.75;
        c_expected(1, 0) = 61.5;
        c_expected(2, 0) = -27.75;
        c_expected(3, 0) = -61.0;
        c_expected(0, 1) = 27.0;
        c_expected(1, 1) = -47.75;
        c_expected(2, 1) = 18.75;
        c_expected(3, 1) = -13.75;
        c_expected(0, 2) = 74.75;
        c_expected(1, 2) = -29.75;
        c_expected(2, 2) = 8.0;
        c_expected(3, 2) = 8.75;
    }

    advdiff::Field c_new{g.nx, g.ny};

    advdiff::apply_periodic_bc(c);
    advdiff::step_seq(c, c_new, g, p, dt);

    for (int j = 0; j < g.ny ; ++j) {
        for (int i = 0; i < g.nx ; ++i) {
            REQUIRE_THAT(c_new(i, j),  Catch::Matchers::WithinRel(c_expected(i, j)));
        }
    }
}

namespace {

double mt19937_output_to_unit(std::uint_fast32_t x) {
    // mt19937 outputs [0, 2^32 - 1]; dividing by 2^32 maps to [0, 1)
    return x / (static_cast<double>(std::mt19937::max()) + 1.0);
}

double uniform(std::mt19937& rng, double lo, double hi) {
    return mt19937_output_to_unit(rng()) * (hi - lo) + lo;
}

} // namespace

TEST_CASE("Applying step_seq matches the naive implementation", "[kernel]") {
    const auto rows = GENERATE(1, 2, 3, 17, 64);
    const auto cols = GENERATE(1, 2, 3, 31);
    const auto rep = GENERATE(range(0, 5));

    const auto seed = Catch::getSeed() + static_cast<std::uint32_t>(rows * 1000003u + cols * 1009u + rep);
    CAPTURE(rows, cols, rep, seed);

    std::mt19937 rng(seed);

    advdiff::Physics p;
    p.u = 2.5;
    p.v = 1.3;
    p.kappa = 0.6;

    advdiff::Grid g;
    g.nx = cols;
    g.ny = rows;

    const auto dt = 0.5;

    advdiff::Field c{g.nx, g.ny};
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            c(i, j) = uniform(rng, -100., 100.);
        }
    }

    advdiff::Field c_new{g.nx, g.ny};

    advdiff::apply_periodic_bc(c);
    advdiff::step_seq(c, c_new, g, p, dt);

    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            double expected = advdiff::details::advdiff_kernel(c(i, j), c(i - 1, j), c(i + 1, j), c(i, j - 1), c(i, j + 1), g, p, dt);
            REQUIRE_THAT(c_new(i, j), Catch::Matchers::WithinRel(expected));
        }
    }
}

TEST_CASE("stability_number is computed properly") {
    advdiff::Physics p;
    advdiff::Grid g;
    double dt = 1e-3;
    double expected = 0.;

    SECTION("u = v = 0, kappa = 0") {
        p.u = 0.;
        p.v = 0.;
        p.kappa = 0.;

        g.nx = GENERATE(5, 10, 100);
        g.ny = GENERATE(5, 10, 100);

        dt = 1e-3;
        expected = 0.;
    }

    SECTION("u = 2, v = 3, kappa = 0.1, dx = 1e-3, dy = 1e-3/1.3") {
        p.u = GENERATE(-2., 2.);
        p.v = GENERATE(-3., 3.);
        p.kappa = .1;

        g.nx = 1e3; // g.nx = g.lx / g.dx
        g.ny = 1.3e3; // g.ny = g.ly / g.dy

        dt = 1e-3;
        expected = 543.9;
    }

    REQUIRE_THAT(advdiff::stability_number(g, p, dt), Catch::Matchers::WithinRel(expected));
}

namespace {

double compute_dt_from_S(const advdiff::Grid& g, const advdiff::Physics& p, double expected_S) {
    return expected_S / advdiff::stability_number(g, p, 1.);
}

} // namespace

TEST_CASE("A uniform field is unchanged", "[kernel]") {
    const auto rows = 64;
    const auto cols = 31;
    const auto c_value = GENERATE(-5., -1e-3, 0., 1e-3, 5.);

    advdiff::Physics p;
    p.u = GENERATE(-13., 0., 5.7);
    p.v = GENERATE(-13., 0., 5.7);
    p.kappa = GENERATE(0., 0.5, 13.);

    advdiff::Grid g;
    g.nx = cols;
    g.ny = rows;

    // u, v, and kappa maximum absolute value is 13
    // In order for S to be less or equal to 1, dt must be sufficiently low.
    const auto max_dt = p.u == 0. and p.v == 0. and p.kappa == 0. ?
                            1e-3 :
                            compute_dt_from_S(g, p, 1.);
    const auto dt = GENERATE_COPY(max_dt/100., max_dt/10., max_dt);

    double c_x = std::abs(p.u) * dt / g.dx();
    double c_y = std::abs(p.v) * dt / g.dy();
    double d_x = p.kappa * dt / (g.dx() * g.dx());
    double d_y = p.kappa * dt / (g.dy() * g.dy());
    double s = c_x + c_y + 2 * d_x + 2 * d_y;

    CAPTURE(rows, cols, c_value, p.u, p.v, p.kappa, dt, s);

    advdiff::Field c{g.nx, g.ny};
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            c(i, j) = c_value;
        }
    }

    advdiff::Field c_new{g.nx, g.ny};

    advdiff::apply_periodic_bc(c);
    advdiff::step_seq(c, c_new, g, p, dt);

    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            REQUIRE(c_new(i, j) == c_value);
        }
    }
}

TEST_CASE("Total mass is conserved", "[kernel]") {
    SKIP("not implemented yet");
}

TEST_CASE("The discrete maximum principle holds when S <= 1", "[kernel]") {
    SKIP("not implemented yet");
}

TEST_CASE("Pure diffusion converges at second order", "[convergence]") {
    SKIP("not implemented yet");
}

TEST_CASE("Pure upwind advection converges at first order", "[convergence]") {
    SKIP("not implemented yet");
}
