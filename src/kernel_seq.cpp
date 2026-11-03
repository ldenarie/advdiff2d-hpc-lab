#include <advdiff/kernel.hpp>
#include <advdiff/details/kernel.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

void advdiff::apply_periodic_bc(Field& c) {
    const int nx = c.nx();
    const int ny = c.ny();
    for (int i = 0; i < nx; ++i) {
        c(i, -1) = c(i, ny - 1);
        c(i, ny) = c(i, 0);
    }
    for (int j = 0; j < ny; ++j) {
        c(-1, j) = c(nx - 1, j);
        c(nx, j) = c(0, j);
    }
}

void advdiff::step_seq(const Field& c, Field& c_new, const Grid& g, const Physics& p, double dt) {
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            c_new(i, j) = details::advdiff_kernel(c(i, j), c(i - 1, j), c(i + 1, j), c(i, j - 1), c(i, j + 1), g, p, dt);
        }
    }
}

double advdiff::stability_number(const Grid& g, const Physics& p, double dt) {
    const double c_x = std::abs(p.u) * dt / g.dx();
    const double c_y = std::abs(p.v) * dt / g.dy();
    const double d_x = p.kappa * dt / (g.dx() * g.dx());
    const double d_y = p.kappa * dt / (g.dy() * g.dy());
    return c_x + c_y + 2 * d_x + 2 * d_y;
}

double advdiff::total_mass(const Field& c, const Grid& g) {
    double sum = 0.;
    for (int j = 0; j < g.ny ; ++j) {
        sum += std::accumulate(c.data() + c.index(0, j), c.data() + c.index(g.nx - 1, j) + 1, 0.);
    }
    return sum * g.dx() * g.dy();
}

double advdiff::field_min(const Field& c) {
    auto min_value = std::numeric_limits<double>::max();
    for (int j = 0; j < c.ny() ; ++j) {
        const auto line_min_element = std::min_element(c.data() + c.index(0, j), c.data() + c.index(c.nx() - 1, j) + 1);
        min_value = std::min(min_value, *line_min_element);
    }
    return min_value;
}

double advdiff::field_max(const Field& c) {
    auto max_value = std::numeric_limits<double>::lowest();
    for (int j = 0; j < c.ny() ; ++j) {
        const auto line_max_element = std::max_element(c.data() + c.index(0, j), c.data() + c.index(c.nx() - 1, j) + 1);
        max_value = std::max(max_value, *line_max_element);
    }
    return max_value;
}