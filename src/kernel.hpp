#pragma once

#include "field.hpp"
#include "problem.hpp"

namespace advdiff {

// Copies the opposite interior edges into the ghost cells (periodic domain).
// Corners are not needed by the 5-point stencil.
void apply_periodic_bc(Field& c);

// One forward Euler step of the upwind / centred scheme (docs/model.md, section 2).
// Reads the interior and ghost cells of `c`; writes the interior of `c_new`.
// Precondition: ghost cells of `c` are up to date.
void step_seq(const Field& c, Field& c_new, const Grid& g, const Physics& p, double dt);

// S = |u| dt/dx + |v| dt/dy + 2 kappa dt/dx^2 + 2 kappa dt/dy^2 (docs/model.md, section 3).
double stability_number(const Grid& g, const Physics& p, double dt);

// Sum of c(i, j) * dx * dy over the interior cells.
double total_mass(const Field& c, const Grid& g);

// Minimum and maximum over the interior cells.
double field_min(const Field& c);
double field_max(const Field& c);

}  // namespace advdiff
