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

void advdiff::step_naive(const Field& c, Field& c_new, const Grid& g, const Physics& p, double dt) {
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            c_new(i, j) = details::advdiff_kernel_naive(c(i, j), c(i - 1, j), c(i + 1, j), c(i, j - 1), c(i, j + 1), g, p, dt);
        }
    }
}

void advdiff::step_weighted(const Field& c, Field& c_new, const Grid& g, const Physics& p, double dt) {
    const details::Constants constants{g, p, dt};
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            c_new(i, j) = constants.kernel_cell_factor * c(i, j) +
                          constants.kernel_west_factor * c(i - 1, j) +
                          constants.kernel_east_factor * c(i + 1, j) +
                          constants.kernel_south_factor * c(i, j - 1) +
                          constants.kernel_north_factor * c(i, j + 1);
        }
    }
}

void advdiff::step_incremental(const Field& c, Field& c_new, const Grid& g, const Physics& p, double dt) {
    const details::Constants constants{g, p, dt};
    for (int j = 0; j < g.ny; ++j) {
        for (int i = 0; i < g.nx; ++i) {
            c_new(i, j) = c(i, j) +
                          constants.kernel_west_factor * (c(i - 1, j) - c(i, j)) +
                          constants.kernel_east_factor * (c(i + 1, j) - c(i, j)) +
                          constants.kernel_south_factor * (c(i, j - 1) - c(i, j)) +
                          constants.kernel_north_factor * (c(i, j + 1) - c(i, j));
        }
    }
}

double advdiff::stability_number(const Grid& g, const Physics& p, double dt) {
    details::Constants constants{g, p, dt};
    return constants.S;
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