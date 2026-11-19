#include <advdiff/kernel.hpp>
#include <advdiff/details/kernel.hpp>

void advdiff::step_incremental_omp(const Field& c, Field& c_new, const Grid& g, const Physics& p, double dt) {
    const details::Constants constants{g, p, dt};
#ifdef _OPENMP
#pragma omp parallel for schedule(static) default(none) shared(g, c, c_new, constants)
#endif
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
