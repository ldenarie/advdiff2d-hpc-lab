#pragma once

#include <advdiff/field.hpp>
#include <advdiff/problem.hpp>

#include <cmath>

namespace advdiff::details {
    inline double advection_term_naive(double cell, double west, double east, double south, double north,
                                 const advdiff::Grid& g, const advdiff::Physics& p) {
        double a_x = p.u < 0 ?
                       -p.u * (east - cell) / g.dx():
                       -p.u * (cell - west) / g.dx();
        double a_y = p.v < 0 ?
                       -p.v * (north - cell) / g.dy():
                       -p.v * (cell - south) / g.dy();
        return a_x + a_y;
    }

    inline double diffusion_term(double cell, double west, double east, double south, double north,
                                 const advdiff::Grid& g, const advdiff::Physics& p) {
        const double x_term = (east - 2 * cell + west) / (g.dx() * g.dx());
        const double y_term = (north - 2 * cell + south) / (g.dy() * g.dy());
        return p.kappa * (x_term + y_term);
    }

    inline double advdiff_kernel_naive(double cell, double west, double east, double south, double north,
                                 const advdiff::Grid& g, const advdiff::Physics& p, double dt) {
        const double a = advection_term_naive(cell, west, east, south, north, g, p);
        const double d = diffusion_term(cell, west, east, south, north, g, p);
        return cell + dt * (a + d);
    }

    struct Constants {
        Constants(const advdiff::Grid& g, const advdiff::Physics& p, double dt) :
                Cx{std::abs(p.u) * dt / g.dx()},
                Cy{std::abs(p.v) * dt / g.dy()},
                Dx{p.kappa * dt / (g.dx() * g.dx())},
                Dy{p.kappa * dt / (g.dy() * g.dy())},
                S{Cx + Cy + 2 * Dx + 2 * Dy},
                kernel_cell_factor{1. - S},
                kernel_west_factor{p.u > 0 ? Cx + Dx : Dx},
                kernel_east_factor{p.u < 0 ? Cx + Dx : Dx},
                kernel_south_factor{p.v > 0 ? Cy + Dy : Dy},
                kernel_north_factor{p.v < 0 ? Cy + Dy : Dy}
        {}

        const double Cx;
        const double Cy;
        const double Dx;
        const double Dy;
        const double S;
        const double kernel_cell_factor;
        const double kernel_west_factor;
        const double kernel_east_factor;
        const double kernel_south_factor;
        const double kernel_north_factor;
    };
}
