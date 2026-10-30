#include <advdiff/kernel.hpp>

#include <stdexcept>

namespace advdiff {

void apply_periodic_bc(Field& c) {
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

// TODO(block 1): implement every function below, following docs/model.md.
// Remove the `(void)` lines and the exceptions as you go.

void step_seq(const Field& c, Field& c_new, const Grid& g, const Physics& p, double dt) {
    (void)c; (void)c_new; (void)g; (void)p; (void)dt;
    throw std::logic_error("step_seq: not implemented yet");
}

double stability_number(const Grid& g, const Physics& p, double dt) {
    (void)g; (void)p; (void)dt;
    throw std::logic_error("stability_number: not implemented yet");
}

double total_mass(const Field& c, const Grid& g) {
    (void)c; (void)g;
    throw std::logic_error("total_mass: not implemented yet");
}

double field_min(const Field& c) {
    (void)c;
    throw std::logic_error("field_min: not implemented yet");
}

double field_max(const Field& c) {
    (void)c;
    throw std::logic_error("field_max: not implemented yet");
}

}  // namespace advdiff
