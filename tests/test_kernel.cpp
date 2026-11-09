// Verification tests. See docs/model.md, section 5, for the expected properties.

#include <advdiff/field.hpp>
#include <advdiff/kernel.hpp>
#include <advdiff/problem.hpp>

#include <advdiff/details/kernel.hpp>

#include <catch2/catch_get_random_seed.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_range.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <numbers>
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

TEST_CASE("Naive kernel function is correct", "[kernel]") {
    using advdiff::details::advdiff_kernel_naive;

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

    REQUIRE_THAT(advdiff_kernel_naive(cell, west, east, south, north, g, p, dt), Catch::Matchers::WithinAbs(expected, 1e-11));
}

TEST_CASE("Kernel factors are correct", "[kernel]") {
    advdiff::Grid g;
    const auto dx = 2.;
    const auto dy = 3.;
    g.nx = 100;
    g.ny = 100;
    g.lx = g.nx * dx;
    g.ly = g.ny * dy;

    const auto dt = 36.;

    advdiff::Physics p;
    double expected_cell_factor;
    double expected_west_factor;
    double expected_east_factor;
    double expected_south_factor;
    double expected_north_factor;

    SECTION("u = -2, v = -3, kappa = 0") {
        p.u = -2.;
        p.v = -3.;
        p.kappa = 0.;

        expected_cell_factor = -71.;
        expected_west_factor = 0.;
        expected_east_factor = 36.;
        expected_north_factor = 36.;
        expected_south_factor = 0.;
    }

    SECTION("u = -2, v = -3, kappa = 1e-3") {
        p.u = -2.;
        p.v = -3.;
        p.kappa = 1e-3;

        expected_cell_factor = -71.026;
        expected_west_factor = 9e-3;
        expected_east_factor = 36.009;
        expected_north_factor = 36.004;
        expected_south_factor = 4e-3;
    }

    SECTION("u = -2, v = -3, kappa = 1") {
        p.u = -2.;
        p.v = -3.;
        p.kappa = 1.;

        expected_cell_factor = -97.;
        expected_west_factor = 9.;
        expected_east_factor = 45.;
        expected_north_factor = 40.;
        expected_south_factor = 4.;
    }

    SECTION("u = -2, v = -3, kappa = 2") {
        p.u = -2.;
        p.v = -3.;
        p.kappa = 2.;

        expected_cell_factor = -123.;
        expected_west_factor = 18.;
        expected_east_factor = 54.;
        expected_north_factor = 44.;
        expected_south_factor = 8.;
    }

    SECTION("u = -2, v = 0, kappa = 0") {
        p.u = -2.;
        p.v = 0.;
        p.kappa = 0.;

        expected_cell_factor = -35.;
        expected_west_factor = 0.;
        expected_east_factor = 36.;
        expected_north_factor = 0.;
        expected_south_factor = 0.;
    }

    SECTION("u = -2, v = 0, kappa = 1e-3") {
        p.u = -2.;
        p.v = 0.;
        p.kappa = 1e-3;

        expected_cell_factor = -35.026;
        expected_west_factor = 9e-3;
        expected_east_factor = 36.009;
        expected_north_factor = 4e-3;
        expected_south_factor = 4e-3;
    }

    SECTION("u = -2, v = 0, kappa = 1") {
        p.u = -2.;
        p.v = 0.;
        p.kappa = 1.;

        expected_cell_factor = -61.;
        expected_west_factor = 9.;
        expected_east_factor = 45.;
        expected_north_factor = 4.;
        expected_south_factor = 4.;
    }

    SECTION("u = -2, v = 0, kappa = 2") {
        p.u = -2.;
        p.v = 0.;
        p.kappa = 2.;

        expected_cell_factor = -87.;
        expected_west_factor = 18.;
        expected_east_factor = 54.;
        expected_north_factor = 8.;
        expected_south_factor = 8.;
    }

    SECTION("u = -2, v = 3, kappa = 0") {
        p.u = -2.;
        p.v = 3.;
        p.kappa = 0.;

        expected_cell_factor = -71.;
        expected_west_factor = 0.;
        expected_east_factor = 36.;
        expected_north_factor = 0.;
        expected_south_factor = 36.;
    }

    SECTION("u = -2, v = 3, kappa = 1e-3") {
        p.u = -2.;
        p.v = 3.;
        p.kappa = 1e-3;

        expected_cell_factor = -71.026;
        expected_west_factor = 9e-3;
        expected_east_factor = 36.009;
        expected_north_factor = 4e-3;
        expected_south_factor = 36.004;
    }

    SECTION("u = -2, v = 3, kappa = 1") {
        p.u = -2.;
        p.v = 3.;
        p.kappa = 1.;

        expected_cell_factor = -97.;
        expected_west_factor = 9.;
        expected_east_factor = 45.;
        expected_north_factor = 4.;
        expected_south_factor = 40.;
    }

    SECTION("u = -2, v = 3, kappa = 2") {
        p.u = -2.;
        p.v = 3.;
        p.kappa = 2.;

        expected_cell_factor = -123.;
        expected_west_factor = 18.;
        expected_east_factor = 54.;
        expected_north_factor = 8.;
        expected_south_factor = 44.;
    }

    SECTION("u = 0, v = -3, kappa = 0") {
        p.u = 0.;
        p.v = -3.;
        p.kappa = 0.;

        expected_cell_factor = -35.;
        expected_west_factor = 0.;
        expected_east_factor = 0.;
        expected_north_factor = 36.;
        expected_south_factor = 0.;
    }

    SECTION("u = 0, v = -3, kappa = 1e-3") {
        p.u = 0.;
        p.v = -3.;
        p.kappa = 1e-3;

        expected_cell_factor = -35.026;
        expected_west_factor = 9e-3;
        expected_east_factor = 9e-3;
        expected_north_factor = 36.004;
        expected_south_factor = 4e-3;
    }

    SECTION("u = 0, v = -3, kappa = 1") {
        p.u = 0.;
        p.v = -3.;
        p.kappa = 1.;

        expected_cell_factor = -61.;
        expected_west_factor = 9.;
        expected_east_factor = 9.;
        expected_north_factor = 40.;
        expected_south_factor = 4.;
    }

    SECTION("u = 0, v = -3, kappa = 2") {
        p.u = 0.;
        p.v = -3.;
        p.kappa = 2.;

        expected_cell_factor = -87.;
        expected_west_factor = 18.;
        expected_east_factor = 18.;
        expected_north_factor = 44.;
        expected_south_factor = 8.;
    }

    SECTION("u = 0, v = 0, kappa = 0") {
        p.u = 0.;
        p.v = 0.;
        p.kappa = 0.;

        expected_cell_factor = 1.;
        expected_west_factor = 0.;
        expected_east_factor = 0.;
        expected_north_factor = 0.;
        expected_south_factor = 0.;
    }

    SECTION("u = 0, v = 0, kappa = 1e-3") {
        p.u = 0.;
        p.v = 0.;
        p.kappa = 1e-3;

        expected_cell_factor = .974;
        expected_west_factor = 9e-3;
        expected_east_factor = 9e-3;
        expected_north_factor = 4e-3;
        expected_south_factor = 4e-3;
    }

    SECTION("u = 0, v = 0, kappa = 1") {
        p.u = 0.;
        p.v = 0.;
        p.kappa = 1.;

        expected_cell_factor = -25.;
        expected_west_factor = 9.;
        expected_east_factor = 9.;
        expected_north_factor = 4.;
        expected_south_factor = 4.;
    }

    SECTION("u = 0, v = 0, kappa = 2") {
        p.u = 0.;
        p.v = 0.;
        p.kappa = 2.;

        expected_cell_factor = -51.;
        expected_west_factor = 18.;
        expected_east_factor = 18.;
        expected_north_factor = 8.;
        expected_south_factor = 8.;
    }

    SECTION("u = 0, v = 3, kappa = 0") {
        p.u = 0.;
        p.v = 3.;
        p.kappa = 0.;

        expected_cell_factor = -35.;
        expected_west_factor = 0.;
        expected_east_factor = 0.;
        expected_north_factor = 0.;
        expected_south_factor = 36.;
    }

    SECTION("u = 0, v = 3, kappa = 1e-3") {
        p.u = 0.;
        p.v = 3.;
        p.kappa = 1e-3;

        expected_cell_factor = -35.026;
        expected_west_factor = 9e-3;
        expected_east_factor = 9e-3;
        expected_north_factor = 4e-3;
        expected_south_factor = 36.004;
    }

    SECTION("u = 0, v = 3, kappa = 1") {
        p.u = 0.;
        p.v = 3.;
        p.kappa = 1.;

        expected_cell_factor = -61.;
        expected_west_factor = 9.;
        expected_east_factor = 9.;
        expected_north_factor = 4.;
        expected_south_factor = 40.;
    }

    SECTION("u = 0, v = 3, kappa = 2") {
        p.u = 0.;
        p.v = 3.;
        p.kappa = 2.;

        expected_cell_factor = -87.;
        expected_west_factor = 18.;
        expected_east_factor = 18.;
        expected_north_factor = 8.;
        expected_south_factor = 44.;
    }

    SECTION("u = 2, v = -3, kappa = 0") {
        p.u = 2.;
        p.v = -3.;
        p.kappa = 0.;

        expected_cell_factor = -71.;
        expected_west_factor = 36.;
        expected_east_factor = 0.;
        expected_north_factor = 36.;
        expected_south_factor = 0.;
    }

    SECTION("u = 2, v = -3, kappa = 1e-3") {
        p.u = 2.;
        p.v = -3.;
        p.kappa = 1e-3;

        expected_cell_factor = -71.026;
        expected_west_factor = 36.009;
        expected_east_factor = 9e-3;
        expected_north_factor = 36.004;
        expected_south_factor = 4e-3;
    }

    SECTION("u = 2, v = -3, kappa = 1") {
        p.u = 2.;
        p.v = -3.;
        p.kappa = 1.;

        expected_cell_factor = -97.;
        expected_west_factor = 45.;
        expected_east_factor = 9.;
        expected_north_factor = 40.;
        expected_south_factor = 4.;
    }

    SECTION("u = 2, v = -3, kappa = 2") {
        p.u = 2.;
        p.v = -3.;
        p.kappa = 2.;

        expected_cell_factor = -123.;
        expected_west_factor = 54.;
        expected_east_factor = 18.;
        expected_north_factor = 44.;
        expected_south_factor = 8.;
    }

    SECTION("u = 2, v = 0, kappa = 0") {
        p.u = 2.;
        p.v = 0.;
        p.kappa = 0.;

        expected_cell_factor = -35.;
        expected_west_factor = 36.;
        expected_east_factor = 0.;
        expected_north_factor = 0.;
        expected_south_factor = 0.;
    }

    SECTION("u = 2, v = 0, kappa = 1e-3") {
        p.u = 2.;
        p.v = 0.;
        p.kappa = 1e-3;

        expected_cell_factor = -35.026;
        expected_west_factor = 36.009;
        expected_east_factor = 9e-3;
        expected_north_factor = 4e-3;
        expected_south_factor = 4e-3;
    }

    SECTION("u = 2, v = 0, kappa = 1") {
        p.u = 2.;
        p.v = 0.;
        p.kappa = 1.;

        expected_cell_factor = -61.;
        expected_west_factor = 45.;
        expected_east_factor = 9.;
        expected_north_factor = 4.;
        expected_south_factor = 4.;
    }

    SECTION("u = 2, v = 0, kappa = 2") {
        p.u = 2.;
        p.v = 0.;
        p.kappa = 2.;

        expected_cell_factor = -87.;
        expected_west_factor = 54.;
        expected_east_factor = 18.;
        expected_north_factor = 8.;
        expected_south_factor = 8.;
    }

    SECTION("u = 2, v = 3, kappa = 0") {
        p.u = 2.;
        p.v = 3.;
        p.kappa = 0.;

        expected_cell_factor = -71.;
        expected_west_factor = 36.;
        expected_east_factor = 0.;
        expected_north_factor = 0.;
        expected_south_factor = 36.;
    }

    SECTION("u = 2, v = 3, kappa = 1e-3") {
        p.u = 2.;
        p.v = 3.;
        p.kappa = 1e-3;

        expected_cell_factor = -71.026;
        expected_west_factor = 36.009;
        expected_east_factor = 9e-3;
        expected_north_factor = 4e-3;
        expected_south_factor = 36.004;
    }

    SECTION("u = 2, v = 3, kappa = 1") {
        p.u = 2.;
        p.v = 3.;
        p.kappa = 1.;

        expected_cell_factor = -97.;
        expected_west_factor = 45.;
        expected_east_factor = 9.;
        expected_north_factor = 4.;
        expected_south_factor = 40.;
    }

    SECTION("u = 2, v = 3, kappa = 2") {
        p.u = 2.;
        p.v = 3.;
        p.kappa = 2.;

        expected_cell_factor = -123.;
        expected_west_factor = 54.;
        expected_east_factor = 18.;
        expected_north_factor = 8.;
        expected_south_factor = 44.;
    }

    CAPTURE(p.u, p.v, p.kappa);

    const advdiff::details::Constants constants{g, p, dt};

    REQUIRE_THAT(expected_cell_factor, Catch::Matchers::WithinRel(expected_cell_factor));
    REQUIRE_THAT(expected_west_factor, Catch::Matchers::WithinRel(expected_west_factor));
    REQUIRE_THAT(expected_east_factor, Catch::Matchers::WithinRel(expected_east_factor));
    REQUIRE_THAT(expected_south_factor, Catch::Matchers::WithinRel(expected_south_factor));
    REQUIRE_THAT(expected_north_factor, Catch::Matchers::WithinRel(expected_north_factor));
}

TEST_CASE("step_naive for specific fields") {
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
    advdiff::step_naive(c, c_new, g, p, dt);

    for (int j = 0; j < g.ny ; ++j) {
        for (int i = 0; i < g.nx ; ++i) {
            REQUIRE_THAT(c_new(i, j),  Catch::Matchers::WithinRel(c_expected(i, j)));
        }
    }
}

TEST_CASE("step_weighted for specific fields") {
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
    advdiff::step_weighted(c, c_new, g, p, dt);

    for (int j = 0; j < g.ny ; ++j) {
        for (int i = 0; i < g.nx ; ++i) {
            REQUIRE_THAT(c_new(i, j),  Catch::Matchers::WithinRel(c_expected(i, j)));
        }
    }
}

TEST_CASE("step_incremental for specific fields") {
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
    advdiff::step_incremental(c, c_new, g, p, dt);

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

TEST_CASE("Applying step_naive matches the naive implementation", "[kernel]") {
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
    advdiff::step_naive(c, c_new, g, p, dt);

    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            double expected = advdiff::details::advdiff_kernel_naive(c(i, j), c(i - 1, j), c(i + 1, j), c(i, j - 1), c(i, j + 1), g, p, dt);
            REQUIRE_THAT(c_new(i, j), Catch::Matchers::WithinRel(expected));
        }
    }
}

TEST_CASE("Applying step_weighted matches the naive implementation", "[kernel]") {
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
    advdiff::step_weighted(c, c_new, g, p, dt);

    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            double expected = advdiff::details::advdiff_kernel_naive(c(i, j), c(i - 1, j), c(i + 1, j), c(i, j - 1), c(i, j + 1), g, p, dt);
            const advdiff::details::Constants k{g, p, dt};
            const double scale = std::abs(k.kernel_cell_factor  * c(i, j))
                   + std::abs(k.kernel_west_factor  * c(i - 1, j))
                   + std::abs(k.kernel_east_factor  * c(i + 1, j))
                   + std::abs(k.kernel_south_factor * c(i, j - 1))
                   + std::abs(k.kernel_north_factor * c(i, j + 1));
            const double tolerance = 16 * std::numeric_limits<double>::epsilon() * scale;
            REQUIRE_THAT(c_new(i, j), Catch::Matchers::WithinAbs(expected, tolerance));
        }
    }
}

TEST_CASE("Applying step_incremental matches the naive implementation", "[kernel]") {
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
    advdiff::step_incremental(c, c_new, g, p, dt);

    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            double expected = advdiff::details::advdiff_kernel_naive(c(i, j), c(i - 1, j), c(i + 1, j), c(i, j - 1), c(i, j + 1), g, p, dt);
            const advdiff::details::Constants k{g, p, dt};
            const double scale = std::abs(k.kernel_cell_factor  * c(i, j))
                   + std::abs(k.kernel_west_factor  * c(i - 1, j))
                   + std::abs(k.kernel_east_factor  * c(i + 1, j))
                   + std::abs(k.kernel_south_factor * c(i, j - 1))
                   + std::abs(k.kernel_north_factor * c(i, j + 1));
            const double tolerance = 16 * std::numeric_limits<double>::epsilon() * scale;
            REQUIRE_THAT(c_new(i, j), Catch::Matchers::WithinAbs(expected, tolerance));
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

TEST_CASE("A uniform field is unchanged with step_weighted", "[kernel]") {
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

    double s = advdiff::stability_number(g, p, dt);

    CAPTURE(rows, cols, c_value, p.u, p.v, p.kappa, dt, s);

    advdiff::Field c{g.nx, g.ny};
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            c(i, j) = c_value;
        }
    }

    advdiff::Field c_new{g.nx, g.ny};

    advdiff::apply_periodic_bc(c);
    advdiff::step_weighted(c, c_new, g, p, dt);

    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            REQUIRE_THAT(c_new(i, j), Catch::Matchers::WithinRel(c_value));
        }
    }
}

TEST_CASE("A uniform field is unchanged with step_incremental", "[kernel]") {
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

    double s = advdiff::stability_number(g, p, dt);

    CAPTURE(rows, cols, c_value, p.u, p.v, p.kappa, dt, s);

    advdiff::Field c{g.nx, g.ny};
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            c(i, j) = c_value;
        }
    }

    advdiff::Field c_new{g.nx, g.ny};

    advdiff::apply_periodic_bc(c);
    advdiff::step_incremental(c, c_new, g, p, dt);

    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            REQUIRE(c_new(i, j) == c_value);
        }
    }
}

TEST_CASE("total_mass properly works") {
    advdiff::Grid g;
    g.nx = 10;
    g.ny = 5;
    advdiff::Field c(g.nx, g.ny);
    for (int j = 0; j < g.ny ; ++j) {
        for (int i = 0; i < g.nx ; ++i) {
            c(i, j) = 1'000 * j + i;
        }
    }
    // Each row sums to 1'000 * nx * j + 45, Total is 1'000 * 100 + 45 * 5 = 100'225
    // Each cell area is 1/50
    REQUIRE_THAT(advdiff::total_mass(c, g), Catch::Matchers::WithinRel(2'004.5));

    // application of the periodic boundary condition does not change the sum
    advdiff::apply_periodic_bc(c);
    REQUIRE_THAT(advdiff::total_mass(c, g), Catch::Matchers::WithinRel(2'004.5));
}

namespace {

struct Gaussian {
    double height;
    double mu_x;
    double mu_y;
    double sigma_x;
    double sigma_y;

    double value(double x, double y) {
        const auto dx = x - mu_x;
        const auto dy = y - mu_y;
        return height * std::exp(-(dx * dx / (2 * sigma_x * sigma_x) + dy * dy / (2 * sigma_y * sigma_y)));
    }
};

} // namespace

TEST_CASE("Total mass is conserved by step_weighted", "[kernel]") {
    const auto rows = 64;
    const auto cols = 31;

    Gaussian gaussian{
        .height = 12.5,
        .mu_x = 0.325,
        .mu_y = 0.87768,
        .sigma_x = 0.1,
        .sigma_y = 0.053
    };

    advdiff::Physics p;
    p.u = GENERATE(-13., 5.7);
    p.v = GENERATE(-13., 5.7);
    p.kappa = GENERATE(0.5, 13.);

    advdiff::Grid g;
    g.nx = cols;
    g.ny = rows;

    const auto dt = compute_dt_from_S(g, p, 1.);

    double s = advdiff::stability_number(g, p, dt);

    CAPTURE(p.u, p.v, p.kappa, dt, s);

    advdiff::Field c{g.nx, g.ny};
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            c(i, j) = gaussian.value(g.x(i), g.y(j));
        }
    }

    advdiff::Field c_new{g.nx, g.ny};

    const auto nb_step = 500;

    const auto expected_total_mass = advdiff::total_mass(c, g);

    for (int step = 0; step < nb_step ; ++step) {
        advdiff::apply_periodic_bc(c);
        advdiff::step_weighted(c, c_new, g, p, dt);
        std::swap(c, c_new);
    }

    REQUIRE_THAT(advdiff::total_mass(c, g), Catch::Matchers::WithinRel(expected_total_mass, 1e-12));
}

TEST_CASE("Total mass is conserved by step_incremental", "[kernel]") {
    const auto rows = 64;
    const auto cols = 31;

    Gaussian gaussian{
        .height = 12.5,
        .mu_x = 0.325,
        .mu_y = 0.87768,
        .sigma_x = 0.1,
        .sigma_y = 0.053
    };

    advdiff::Physics p;
    p.u = GENERATE(-13., 5.7);
    p.v = GENERATE(-13., 5.7);
    p.kappa = GENERATE(0.5, 13.);

    advdiff::Grid g;
    g.nx = cols;
    g.ny = rows;

    const auto dt = compute_dt_from_S(g, p, 1.);

    double s = advdiff::stability_number(g, p, dt);

    CAPTURE(p.u, p.v, p.kappa, dt, s);

    advdiff::Field c{g.nx, g.ny};
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            c(i, j) = gaussian.value(g.x(i), g.y(j));
        }
    }

    advdiff::Field c_new{g.nx, g.ny};

    const auto nb_step = 500;

    const auto expected_total_mass = advdiff::total_mass(c, g);

    for (int step = 0; step < nb_step ; ++step) {
        advdiff::apply_periodic_bc(c);
        advdiff::step_incremental(c, c_new, g, p, dt);
        std::swap(c, c_new);
    }

    REQUIRE_THAT(advdiff::total_mass(c, g), Catch::Matchers::WithinRel(expected_total_mass, 1e-12));
}

TEST_CASE("field_min and field_max works") {
    advdiff::Grid g;
    g.nx = 4;
    g.ny = 3;
    advdiff::Field c(g.nx, g.ny);

    double expected_min;
    double expected_max;

    SECTION ("All strictly positive") {
        c(0, 0) = 1.;
        c(1, 0) = 5.;
        c(2, 0) = 13.;
        c(3, 0) = 17.;
        c(0, 1) = 12.;
        c(1, 1) = 8.;
        c(2, 1) = 3.;
        c(3, 1) = 4.;
        c(0, 2) = 19.;
        c(1, 2) = 2.;
        c(2, 2) = 2.;
        c(3, 2) = 4.;

        expected_min = 1.;
        expected_max = 19.;
    }

    SECTION ("All strictly negative") {
        c(0, 0) = -1.;
        c(1, 0) = -5.;
        c(2, 0) = -13.;
        c(3, 0) = -17.;
        c(0, 1) = -12.;
        c(1, 1) = -8.;
        c(2, 1) = -3.;
        c(3, 1) = -4.;
        c(0, 2) = -19.;
        c(1, 2) = -2.;
        c(2, 2) = -2.;
        c(3, 2) = -4.;

        expected_min = -19.;
        expected_max = -1.;
    }

    SECTION ("All positive") {
        c(0, 0) = 1.;
        c(1, 0) = 5.;
        c(2, 0) = 13.;
        c(3, 0) = 17.;
        c(0, 1) = 12.;
        c(1, 1) = 8.;
        c(2, 1) = 3.;
        c(3, 1) = 4.;
        c(0, 2) = 19.;
        c(1, 2) = 0.;
        c(2, 2) = 2.;
        c(3, 2) = 4.;

        expected_min = 0.;
        expected_max = 19.;
    }

    SECTION ("All negative") {
        c(0, 0) = -1.;
        c(1, 0) = -5.;
        c(2, 0) = -13.;
        c(3, 0) = -17.;
        c(0, 1) = -12.;
        c(1, 1) = -8.;
        c(2, 1) = -3.;
        c(3, 1) = -4.;
        c(0, 2) = -19.;
        c(1, 2) = 0.;
        c(2, 2) = -2.;
        c(3, 2) = -4.;

        expected_min = -19.;
        expected_max = 0.;
    }

    SECTION ("Negative and positive") {
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

        expected_min = -19.;
        expected_max = 17.;
    }

    REQUIRE_THAT(advdiff::field_min(c), Catch::Matchers::WithinRel(expected_min));
    REQUIRE_THAT(advdiff::field_max(c), Catch::Matchers::WithinRel(expected_max));

    // application of the periodic boundary condition does not change the min and max
    advdiff::apply_periodic_bc(c);
    REQUIRE_THAT(advdiff::field_min(c), Catch::Matchers::WithinRel(expected_min));
    REQUIRE_THAT(advdiff::field_max(c), Catch::Matchers::WithinRel(expected_max));
}

TEST_CASE("The discrete maximum principle holds when S <= 1 with step_weighted", "[kernel]") {
    const auto rows = 64;
    const auto cols = 31;

    Gaussian gaussian{
        .height = 12.5,
        .mu_x = 0.325,
        .mu_y = 0.87768,
        .sigma_x = 0.1,
        .sigma_y = 0.053
    };

    advdiff::Physics p;
    p.u = GENERATE(-13., 5.7);
    p.v = GENERATE(-13., 5.7);
    p.kappa = GENERATE(0.5, 13.);

    advdiff::Grid g;
    g.nx = cols;
    g.ny = rows;

    const auto dt = compute_dt_from_S(g, p, 1.);

    double s = advdiff::stability_number(g, p, dt);

    CAPTURE(p.u, p.v, p.kappa, dt, s);

    advdiff::Field c{g.nx, g.ny};
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            c(i, j) = gaussian.value(g.x(i), g.y(j));
        }
    }

    advdiff::Field c_new{g.nx, g.ny};

    const auto nb_step = 500;

    const auto initial_min_value = advdiff::field_min(c);
    const auto initial_max_value = advdiff::field_max(c);
    REQUIRE(initial_min_value <= initial_max_value);

    for (int step = 0; step < nb_step ; ++step) {
        advdiff::apply_periodic_bc(c);
        advdiff::step_weighted(c, c_new, g, p, dt);
        std::swap(c, c_new);
    }

    const auto final_min_value = advdiff::field_min(c);
    const auto final_max_value = advdiff::field_max(c);

    REQUIRE(final_min_value <= final_max_value);
    REQUIRE(final_min_value >= initial_min_value - 1e-14);
    REQUIRE(final_max_value <= initial_max_value + 1e-14);
}

TEST_CASE("The discrete maximum principle holds when S <= 1 with step_incremental", "[kernel]") {
    const auto rows = 64;
    const auto cols = 31;

    Gaussian gaussian{
        .height = 12.5,
        .mu_x = 0.325,
        .mu_y = 0.87768,
        .sigma_x = 0.1,
        .sigma_y = 0.053
    };

    advdiff::Physics p;
    p.u = GENERATE(-13., 5.7);
    p.v = GENERATE(-13., 5.7);
    p.kappa = GENERATE(0.5, 13.);

    advdiff::Grid g;
    g.nx = cols;
    g.ny = rows;

    const auto dt = compute_dt_from_S(g, p, 1.);

    double s = advdiff::stability_number(g, p, dt);

    CAPTURE(p.u, p.v, p.kappa, dt, s);

    advdiff::Field c{g.nx, g.ny};
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            c(i, j) = gaussian.value(g.x(i), g.y(j));
        }
    }

    advdiff::Field c_new{g.nx, g.ny};

    const auto nb_step = 500;

    const auto initial_min_value = advdiff::field_min(c);
    const auto initial_max_value = advdiff::field_max(c);
    REQUIRE(initial_min_value <= initial_max_value);

    for (int step = 0; step < nb_step ; ++step) {
        advdiff::apply_periodic_bc(c);
        advdiff::step_incremental(c, c_new, g, p, dt);
        std::swap(c, c_new);
    }

    const auto final_min_value = advdiff::field_min(c);
    const auto final_max_value = advdiff::field_max(c);

    REQUIRE(final_min_value <= final_max_value);
    REQUIRE(final_min_value >= initial_min_value - 1e-14);
    REQUIRE(final_max_value <= initial_max_value + 1e-14);
}

namespace {

template <class ExactSolutionFunction>
double rms_error(const advdiff::Field& c, const advdiff::Grid& g, double T, ExactSolutionFunction&& exact_solution) {
    double total_square_error = 0.;
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            const auto cell_error = c(i, j) - exact_solution(g.x(i), g.y(j), T);
            total_square_error += (cell_error * cell_error);
        }
    }
    return std::sqrt(total_square_error / (g.nx * g.ny));
}

template <class ComputeDeltaT, class ExactSolutionFunction>
double convergence_rms_error(std::string_view kernel, int n, double u, double v, double kappa, double T,
                             ComputeDeltaT&& dt_computer, ExactSolutionFunction&& exact_solution) {
    auto kernel_function = [kernel](){
        if (kernel == "weighted") {
            return advdiff::step_weighted;
        } else {
            return advdiff::step_incremental;
        }
    }();
    advdiff::Physics p{
        .u = u,
        .v = v,
        .kappa = kappa
    };

    advdiff::Grid g {
        .nx = n,
        .ny = n
    };

    advdiff::Field c{g.nx, g.ny};
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            c(i, j) = exact_solution(g.x(i), g.y(j), 0.);
        }
    }

    advdiff::Field c_new{g.nx, g.ny};

    const auto nb_step = static_cast<int>(std::ceil(T / dt_computer(g, p)));
    const auto dt = T / nb_step;

    for (int step = 0; step < nb_step ; ++step) {
        advdiff::apply_periodic_bc(c);
        kernel_function(c, c_new, g, p, dt);
        std::swap(c, c_new);
    }

    return rms_error(c, g, T, exact_solution);
}

} // namespace

TEST_CASE("Pure diffusion converges at second order", "[convergence]") {
    const auto kappa = .1;
    const auto T = .05;

    auto pure_diffusion_exact_solution = [kappa](double x, double y, double t) {
        using std::numbers::pi;
        return 1 +
                  std::sin(2 * pi * x) *
                  std::sin(2 * pi * y) *
                  std::exp(-8 * pi * pi * kappa * t);
    };

    auto dt_computer = [](const advdiff::Grid& g, const advdiff::Physics& p) {
        return 0.2 * g.dx() * g.dx() / p.kappa;
    };

    auto kernel = GENERATE("weighted", "incremental");

    const auto e_32 = convergence_rms_error(kernel, 32, 0., 0., kappa, T, dt_computer, pure_diffusion_exact_solution);
    const auto e_64 = convergence_rms_error(kernel, 64, 0., 0., kappa, T, dt_computer, pure_diffusion_exact_solution);
    const auto e_128 = convergence_rms_error(kernel, 128, 0., 0., kappa, T, dt_computer, pure_diffusion_exact_solution);

    const auto order_32_64 = std::log2(e_32 / e_64);
    const auto order_64_128 = std::log2(e_64 / e_128);
    CAPTURE(kernel, e_32, e_64, e_128, order_32_64, order_64_128);

    REQUIRE_THAT(order_32_64, Catch::Matchers::WithinAbs(2., .2));
    REQUIRE_THAT(order_64_128, Catch::Matchers::WithinAbs(2., .2));
}

TEST_CASE("Pure upwind advection converges at first order", "[convergence]") {
    const auto u = 1.;
    const auto v = .5;
    const auto T = .25;

    auto pure_advection_exact_solution = [u, v](double x, double y, double t) {
        using std::numbers::pi;
        return 1 +
                  std::sin(2 * pi * (x - u * t)) *
                  std::sin(2 * pi * (y - v * t));
    };

    auto dt_computer = [](const advdiff::Grid& g, const advdiff::Physics& p) {
        return compute_dt_from_S(g, p, .5);
    };

    auto kernel = GENERATE("weighted", "incremental");

    const auto e_64 = convergence_rms_error(kernel, 64, u, v, 0., T, dt_computer, pure_advection_exact_solution);
    const auto e_128 = convergence_rms_error(kernel, 128, u, v, 0., T, dt_computer, pure_advection_exact_solution);
    const auto e_256 = convergence_rms_error(kernel, 256, u, v, 0., T, dt_computer, pure_advection_exact_solution);

    const auto order_64_128 = std::log2(e_64 / e_128);
    const auto order_128_256 = std::log2(e_128 / e_256);
    CAPTURE(kernel, e_64, e_128, e_256, order_64_128, order_128_256);

    REQUIRE_THAT(order_64_128, Catch::Matchers::WithinAbs(1., .1));
    REQUIRE_THAT(order_128_256, Catch::Matchers::WithinAbs(1., .1));
}
