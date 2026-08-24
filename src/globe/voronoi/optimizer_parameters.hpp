#ifndef GLOBEART_SRC_GLOBE_VORONOI_OPTIMIZER_PARAMETERS_HPP_
#define GLOBEART_SRC_GLOBE_VORONOI_OPTIMIZER_PARAMETERS_HPP_

#include <cstddef>
#include <string>

namespace globe::voronoi {

struct NewtonParameters {
    size_t max_iterations = 100;
    // The gradient is a difference of moment integrals, so it bottoms out
    // near the square root of machine epsilon like the capacity error does.
    double gradient_tolerance = 1e-7;
    double initial_trust_radius = 0.1;
    double max_trust_radius = 1.0;
    double minimum_trust_radius = 1e-12;
    double acceptance_threshold = 0.1;
    size_t max_conjugate_gradient_iterations = 100;
    double conjugate_gradient_tolerance = 0.1;
    std::string curvature = "exact";
};

struct NewtonReport {
    size_t iterations = 0;
    size_t accepted_steps = 0;
    double cvt_energy = 0.0;
    double gradient_norm = 0.0;
    bool converged = false;
    bool stalled = false;
};

struct CapacityConstrainedParameters {
    size_t max_outer_iterations = 30;
    size_t max_inner_iterations = 200;
    double relative_capacity_tolerance = 1e-7;
    double penalty_growth = 10.0;
    double required_violation_decrease = 0.25;
    double max_penalty_growth_factor = 1e8;
    size_t max_stalled_outer_iterations = 2;
    size_t lbfgs_history = 8;
    std::string inner_solver = "lbfgs";
    NewtonParameters newton;
};

struct CapacityConstrainedReport {
    size_t outer_iterations = 0;
    size_t inner_iterations = 0;
    double cvt_energy = 0.0;
    double relative_rms_capacity_error = 0.0;
    bool converged = false;
    bool stalled = false;
};

} // namespace globe::voronoi

#endif //GLOBEART_SRC_GLOBE_VORONOI_OPTIMIZER_PARAMETERS_HPP_
