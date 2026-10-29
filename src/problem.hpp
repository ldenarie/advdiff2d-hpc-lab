#pragma once

namespace advdiff {

// Uniform periodic grid on [0, lx) x [0, ly), cell-centred values.
struct Grid {
    int nx;
    int ny;
    double lx = 1.0;
    double ly = 1.0;

    double dx() const { return lx / nx; }
    double dy() const { return ly / ny; }
    double x(int i) const { return (i + 0.5) * dx(); }
    double y(int j) const { return (j + 0.5) * dy(); }
};

// Uniform velocity (u, v) and diffusion coefficient kappa.
struct Physics {
    double u = 0.0;
    double v = 0.0;
    double kappa = 0.0;
};

}  // namespace advdiff
