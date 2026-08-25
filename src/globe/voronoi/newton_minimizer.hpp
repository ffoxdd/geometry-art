#ifndef GLOBEART_SRC_GLOBE_VORONOI_NEWTON_MINIMIZER_HPP_
#define GLOBEART_SRC_GLOBE_VORONOI_NEWTON_MINIMIZER_HPP_

#include "optimizer_parameters.hpp"
#include "trust_region_step.hpp"
#include "../types.hpp"
#include <chrono>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <optional>
#include <utility>
#include <vector>

namespace globe::voronoi {

// What a geometry must supply for the trust-region Newton loop to descend
// its Lagrangian: the value and gradient at the current iterate, a
// curvature operator there, and the ability to evaluate a trial step and
// commit to it. The gradient and the step live in whatever parameterisation
// the geometry descends -- tangent planes on the sphere, the plane itself
// on flat domains.
template<typename T>
concept NewtonModel = requires(T& model, const std::vector<Vector3>& step, typename T::Trial&& trial) {
    { model.value() } -> std::convertible_to<double>;
    { model.gradient() } -> std::convertible_to<std::vector<Vector3>>;
    model.refresh_curvature();
    { model.curvature() };
    { model.trial(step) } -> std::convertible_to<typename T::Trial>;
    { model.trial_value(std::declval<const typename T::Trial&>()) } -> std::convertible_to<double>;
    model.accept(std::move(trial));
};

// A long inner solve reports nothing until its outer iteration ends, which
// at large site counts is many minutes of silence, so the descent prints on
// a wall-clock cadence. A solve that finishes sooner stays silent.
constexpr std::chrono::seconds NEWTON_PROGRESS_INTERVAL{5};

// One trust-region Newton descent, shared by every geometry: the loop owns
// the radius policy and the accept/reject protocol, the model owns all the
// geometry. A rejected step leaves the iterate where it was, so the
// curvature there is reused rather than reassembled.
template<NewtonModel ModelType>
size_t newton_minimize(ModelType& model, const NewtonParameters& parameters, size_t max_iterations) {
    TrustRegionStep solver(
        parameters.max_conjugate_gradient_iterations,
        parameters.conjugate_gradient_tolerance
    );

    double radius = parameters.initial_trust_radius;
    size_t iterations = 0;
    bool curvature_stale = true;
    auto last_report = std::chrono::steady_clock::now();

    while (iterations < max_iterations) {
        std::vector<Vector3> gradient = model.gradient();
        double gradient_norm = 0.0;

        for (const Vector3& entry : gradient) {
            gradient_norm += entry.squaredNorm();
        }

        if (std::sqrt(gradient_norm) <= parameters.gradient_tolerance) {
            break;
        }

        auto now = std::chrono::steady_clock::now();

        if (now - last_report >= NEWTON_PROGRESS_INTERVAL) {
            last_report = now;

            std::cout << "    " << std::setw(6) << std::left << "newton" << std::right <<
                std::setw(4) << iterations << ": value " <<
                std::scientific << std::setprecision(6) << model.value() <<
                ", gradient " << std::setprecision(3) << std::sqrt(gradient_norm) <<
                std::defaultfloat << std::endl;
        }

        if (curvature_stale) {
            model.refresh_curvature();
            curvature_stale = false;
        }

        TrustRegionStep::Result step = solver.solve(gradient, model.curvature(), radius);
        ++iterations;

        if (step.predicted_decrease <= 0.0) {
            break;
        }

        typename ModelType::Trial trial = model.trial(step.step);
        double ratio = (model.value() - model.trial_value(trial)) / step.predicted_decrease;

        if (ratio < 0.25) {
            radius *= 0.25;
        } else if (ratio > 0.75 && step.hit_boundary) {
            radius = std::min(2.0 * radius, parameters.max_trust_radius);
        }

        if (ratio <= parameters.acceptance_threshold) {
            if (radius < parameters.minimum_trust_radius) {
                break;
            }

            continue;
        }

        model.accept(std::move(trial));
        curvature_stale = true;
    }

    return iterations;
}

} // namespace globe::voronoi

#endif //GLOBEART_SRC_GLOBE_VORONOI_NEWTON_MINIMIZER_HPP_
