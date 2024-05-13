#pragma once
#include <cmath>
#include "point.hpp"
#include <algorithm>

#define M_UTILS_EPSILON 1e-5
#define M_BOUNDARY (1.0 - M_UTILS_EPSILON)
#define M_MACHINE_EPS std::numeric_limits<double>::epsilon()

namespace hyperbolic_utils {

    inline double sq_norm(const double& a, const double& b) {
        return a * a + b * b;
    }

    inline double poincare_to_klein(const double& a, const double& sq_n) {
        return 2 * a / (1 + sq_n);
    }

    inline double klein_to_poincare(const double& a, const double& sq_n) {
        return a / (1 + std::sqrt(1 - sq_n));
    }

    inline double lorentz_factor(const double& sq_n) {
        return 1 / std::sqrt(1 - sq_n);
    }

    static bool isBoxWithinUnitCircle(const double& min_bounds_x, const double& min_bounds_y, const double& max_bounds_x, const double& max_bounds_y) {
        double d1 = min_bounds_x*min_bounds_x;
        double d2 = min_bounds_y*min_bounds_y;
        if (d1 + d2 < 1.0) {
            double d3 = max_bounds_x*max_bounds_x;
            if (d2 + d3 < 1.0) {
                double d4 = max_bounds_y*max_bounds_y;
                if(d1 + d4 < 1.0 && d3 + d4 < 1.0)
                    return true;
            }
        }
        return false;
    }

    void distance_grad(double u0, double u1, double v0, double v1, double& res_grad_x, double& res_grad_y) {
        if (fabs(u0 - v0) <= M_UTILS_EPSILON && fabs(u1 - v1) <= M_UTILS_EPSILON) {
            res_grad_x = 0;
            res_grad_x = 0;
        }

        double a = u0 - v0;
        double b = u1 - v1;
        double uv2 = a * a + b * b;

        double u_sq = std::clamp(u0 * u0 + u1 * u1, 0.0, M_BOUNDARY);
        double v_sq = std::clamp(v0 * v0 + v1 * v1, 0.0, M_BOUNDARY);
        double alpha = 1. - u_sq;
        double beta = 1. - v_sq;

        double gamma = 1. + (2. / (alpha * beta)) * uv2;
        double shared_scalar = 4. / std::fmax(beta * sqrt((gamma * gamma) - 1.), M_MACHINE_EPS);

        double u_scalar = (v_sq - 2. * (u0 * v0 + u1 * v1) + 1.) / (alpha * alpha);
        double v_scalar = 1. / alpha;

        res_grad_x = shared_scalar * (u_scalar * u0 - v_scalar * v0);
        res_grad_y = shared_scalar * (u_scalar * u1 - v_scalar * v1);
    }
}