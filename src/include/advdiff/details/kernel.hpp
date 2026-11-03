#pragma once

#include <advdiff/field.hpp>
#include <advdiff/problem.hpp>

namespace advdiff::details {
    inline double advection_term(double cell, double west, double east, double south, double north,
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

    inline double advdiff_kernel(double cell, double west, double east, double south, double north,
                                 const advdiff::Grid& g, const advdiff::Physics& p, double dt) {
        const double a = advection_term(cell, west, east, south, north, g, p);
        const double d = diffusion_term(cell, west, east, south, north, g, p);
        return cell + dt * (a + d);
    }
}