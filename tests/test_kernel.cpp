// Verification tests. See docs/model.md, section 5, for the expected properties.
//
// TODO(block 2): implement each test and remove its SKIP().

#include <advdiff/field.hpp>
#include <advdiff/kernel.hpp>
#include <advdiff/problem.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

namespace {

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

TEST_CASE("A uniform field is unchanged", "[kernel]") {
    SKIP("not implemented yet");
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

}  // namespace
