// Verification tests. See docs/model.md, section 5, for the expected properties.
//
// TODO(block 2): implement each test and remove its SKIP().

#include <advdiff/field.hpp>
#include <advdiff/kernel.hpp>
#include <advdiff/problem.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

namespace {

using advdiff::Field;
using advdiff::Grid;
using advdiff::Physics;

TEST_CASE("Ghost cells copy the opposite interior edges", "[bc]") {
    SKIP("not implemented yet");
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
