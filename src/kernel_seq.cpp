#include "kernel.hpp"

#include <stdexcept>

namespace advdiff {

// TODO(block 1): implement every function below, following docs/model.md.
// Remove the `(void)` lines and the exceptions as you go.

void apply_periodic_bc(Field& c) {
    (void)c;
    throw std::logic_error("apply_periodic_bc: not implemented yet");
}

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
