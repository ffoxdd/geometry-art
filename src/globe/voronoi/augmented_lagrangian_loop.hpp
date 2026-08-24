#ifndef GLOBEART_SRC_GLOBE_VORONOI_AUGMENTED_LAGRANGIAN_LOOP_HPP_
#define GLOBEART_SRC_GLOBE_VORONOI_AUGMENTED_LAGRANGIAN_LOOP_HPP_

#include "lagrangian_evaluation.hpp"
#include "optimizer_parameters.hpp"
#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace globe::voronoi {

// What a geometry must supply for the augmented Lagrangian outer loop: an
// inner minimisation at the current multipliers and penalty, and an
// evaluation of the Lagrangian there.
template<typename T>
concept ConstrainedProblem = requires(
    T& problem,
    const std::vector<double>& multipliers,
    double penalty
) {
    { problem.minimize(multipliers, penalty) } -> std::convertible_to<size_t>;
    { problem.evaluate(multipliers, penalty) } -> std::convertible_to<LagrangianEvaluation>;
    { problem.target_mass() } -> std::convertible_to<double>;
    { problem.site_count() } -> std::convertible_to<size_t>;
};

// The outer loop, shared by every geometry: it owns the penalty schedule,
// the multiplier updates and the stall detection, and reports the same way
// everywhere. The problem owns all the geometry.
template<ConstrainedProblem ProblemType>
CapacityConstrainedReport augmented_lagrangian_loop(
    ProblemType& problem,
    const CapacityConstrainedParameters& parameters,
    const std::string& inner_name
) {
    CapacityConstrainedReport report;
    std::vector<double> multipliers(problem.site_count(), 0.0);

    LagrangianEvaluation initial = problem.evaluate(multipliers, 0.0);
    double squared_violation = 0.0;

    for (double error : initial.capacity_errors) {
        squared_violation += error * error;
    }

    double target = problem.target_mass();
    double floor = static_cast<double>(problem.site_count()) *
        std::pow(parameters.relative_capacity_tolerance * target, 2);
    double initial_penalty = 2.0 * std::max(initial.cvt_energy, std::numeric_limits<double>::min()) /
        std::max(squared_violation, floor);
    double penalty = initial_penalty;

    double previous_violation = std::numeric_limits<double>::infinity();
    size_t stalled_iterations = 0;

    for (size_t outer = 0; outer < parameters.max_outer_iterations; ++outer) {
        size_t inner_iterations = problem.minimize(multipliers, penalty);
        LagrangianEvaluation evaluation = problem.evaluate(multipliers, penalty);

        report.outer_iterations = outer + 1;
        report.inner_iterations += inner_iterations;
        report.cvt_energy = evaluation.cvt_energy;
        report.relative_rms_capacity_error = evaluation.root_mean_square_capacity_error() / target;

        std::cout << "  " << std::setw(8) << std::left << "CCVT" << std::right <<
            std::setw(3) << outer + 1 << " (" << std::setw(4) << inner_iterations << " " << inner_name << ")" <<
            ": CVT energy " << std::scientific << std::setprecision(6) << evaluation.cvt_energy <<
            ", capacity RMS " << std::setprecision(3) << report.relative_rms_capacity_error <<
            ", penalty " << std::setprecision(2) << penalty <<
            std::defaultfloat << std::endl;

        if (report.relative_rms_capacity_error < parameters.relative_capacity_tolerance) {
            report.converged = true;
            break;
        }

        double violation = evaluation.max_absolute_capacity_error();
        stalled_iterations = violation < previous_violation ? 0 : stalled_iterations + 1;

        if (stalled_iterations >= parameters.max_stalled_outer_iterations) {
            report.stalled = true;
            break;
        }

        for (size_t i = 0; i < multipliers.size(); ++i) {
            multipliers[i] += penalty * evaluation.capacity_errors[i];
        }

        bool insufficient_decrease = violation > parameters.required_violation_decrease * previous_violation;
        bool below_cap = penalty * parameters.penalty_growth <=
            initial_penalty * parameters.max_penalty_growth_factor;

        if (insufficient_decrease && below_cap) {
            penalty *= parameters.penalty_growth;
        }

        previous_violation = violation;
    }

    return report;
}

} // namespace globe::voronoi

#endif //GLOBEART_SRC_GLOBE_VORONOI_AUGMENTED_LAGRANGIAN_LOOP_HPP_
