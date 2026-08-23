#ifndef GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_CAPACITY_CONSTRAINED_OPTIMIZER_HPP_
#define GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_CAPACITY_CONSTRAINED_OPTIMIZER_HPP_

#include "capacity_constrained_lagrangian.hpp"
#include "capacity_constrained_hessian.hpp"
#include "cvt_hessian.hpp"
#include "newton_optimizer/finite_difference_hessian.hpp"
#include "newton_optimizer/newton_optimizer.hpp"
#include "../../../math/normalization.hpp"
#include "../../../types.hpp"
#include "../../../fields/spherical/field.hpp"
#include "../../../fields/spherical/polynomial_field.hpp"
#include "../core/sphere.hpp"
#include "../core/callback.hpp"
#include <Eigen/Core>
#include <LBFGS.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace globe::voronoi::spherical {

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

template<fields::spherical::Field FieldType = fields::spherical::PolynomialField>
class CapacityConstrainedOptimizer {
 public:
    CapacityConstrainedOptimizer(
        std::unique_ptr<Sphere> sphere,
        FieldType field,
        CapacityConstrainedParameters parameters,
        Callback callback
    );

    std::unique_ptr<Sphere> optimize();
    [[nodiscard]] const CapacityConstrainedReport& report() const { return _report; }

 private:
    class Objective {
     public:
        Objective(CapacityConstrainedOptimizer& optimizer);
        double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& gradient);

        [[nodiscard]] const Eigen::VectorXd& best_x() const { return _best_x; }
        [[nodiscard]] size_t evaluations() const { return _evaluations; }

     private:
        CapacityConstrainedOptimizer& _optimizer;
        Eigen::VectorXd _best_x;
        double _best_value = std::numeric_limits<double>::infinity();
        size_t _evaluations = 0;
    };

    std::unique_ptr<Sphere> _sphere;
    CapacityConstrainedLagrangian<FieldType> _lagrangian;
    CvtHessian<FieldType> _curvature;
    CapacityConstrainedParameters _parameters;
    Callback _callback;
    std::vector<double> _multipliers;
    double _penalty = 0.0;
    double _initial_penalty = 0.0;
    CapacityConstrainedReport _report;

    static constexpr double FINITE_DIFFERENCE_DISPLACEMENT = 1e-5;

    [[nodiscard]] double initial_penalty() const;
    size_t minimize_lagrangian();
    [[nodiscard]] LagrangianEvaluation evaluate() const;
    void update_multipliers(const LagrangianEvaluation& evaluation, double previous_violation);
    void print_progress(size_t outer_iteration, size_t inner_iterations, const LagrangianEvaluation& evaluation) const;

    [[nodiscard]] size_t minimize_lagrangian_by_lbfgs();
    [[nodiscard]] size_t minimize_lagrangian_by_newton();

    struct CurvatureOperator {
        std::function<std::vector<Vector3>(const std::vector<Vector3>&)> apply;
        [[nodiscard]] std::vector<Vector3> multiply(const std::vector<Vector3>& directions) const { return apply(directions); }
    };

    [[nodiscard]] CurvatureOperator curvature_operator(
        const std::vector<Vector3>& points,
        const SphereState& state,
        const LagrangianEvaluation& evaluation,
        const std::vector<Vector3>& gradient
    ) const;
    [[nodiscard]] std::vector<Vector3> tangential_gradient_at(const std::vector<Vector3>& points) const;
    [[nodiscard]] std::vector<Vector3> site_points() const;
    void apply_points(const std::vector<Vector3>& points);

    [[nodiscard]] static std::vector<Vector3> tangential_gradient(
        const std::vector<Vector3>& site_gradients,
        const std::vector<Vector3>& points
    );
    [[nodiscard]] static double gradient_norm(const std::vector<Vector3>& gradient);
    [[nodiscard]] static std::vector<Vector3> stepped(
        const std::vector<Vector3>& points,
        const std::vector<Vector3>& step
    );
    [[nodiscard]] Eigen::VectorXd sites_to_vector() const;
    void apply_sites(const Eigen::VectorXd& x);
    [[nodiscard]] Eigen::VectorXd chain_rule_gradient(const Eigen::VectorXd& x, const std::vector<Vector3>& site_gradients) const;
};

template<fields::spherical::Field FieldType>
CapacityConstrainedOptimizer<FieldType>::CapacityConstrainedOptimizer(
    std::unique_ptr<Sphere> sphere,
    FieldType field,
    CapacityConstrainedParameters parameters,
    Callback callback
) :
    _sphere(std::move(sphere)),
    _lagrangian(field, field.total_mass() / static_cast<double>(_sphere->size())),
    _curvature(field),
    _parameters(parameters),
    _callback(std::move(callback)),
    _multipliers(_sphere->size(), 0.0) {
}

template<fields::spherical::Field FieldType>
std::unique_ptr<Sphere> CapacityConstrainedOptimizer<FieldType>::optimize() {
    _initial_penalty = initial_penalty();
    _penalty = _initial_penalty;
    double previous_violation = std::numeric_limits<double>::infinity();
    size_t stalled_iterations = 0;
    _report = CapacityConstrainedReport{};

    for (size_t outer = 0; outer < _parameters.max_outer_iterations; ++outer) {
        size_t inner_iterations = minimize_lagrangian();
        LagrangianEvaluation evaluation = evaluate();

        _report.outer_iterations = outer + 1;
        _report.inner_iterations += inner_iterations;
        _report.cvt_energy = evaluation.cvt_energy;
        _report.relative_rms_capacity_error = evaluation.root_mean_square_capacity_error() / _lagrangian.target_mass();
        print_progress(outer + 1, inner_iterations, evaluation);

        if (_report.relative_rms_capacity_error < _parameters.relative_capacity_tolerance) {
            _report.converged = true;
            break;
        }

        double violation = evaluation.max_absolute_capacity_error();
        stalled_iterations = violation < previous_violation ? 0 : stalled_iterations + 1;

        if (stalled_iterations >= _parameters.max_stalled_outer_iterations) {
            _report.stalled = true;
            break;
        }

        update_multipliers(evaluation, previous_violation);
        previous_violation = violation;
    }

    return std::move(_sphere);
}

template<fields::spherical::Field FieldType>
double CapacityConstrainedOptimizer<FieldType>::initial_penalty() const {
    LagrangianEvaluation evaluation = evaluate();
    double squared_violation = 0.0;

    for (double error : evaluation.capacity_errors) {
        squared_violation += error * error;
    }

    double target = _lagrangian.target_mass();
    double floor = static_cast<double>(_sphere->size()) * std::pow(_parameters.relative_capacity_tolerance * target, 2);
    return 2.0 * std::max(evaluation.cvt_energy, std::numeric_limits<double>::min()) / std::max(squared_violation, floor);
}

template<fields::spherical::Field FieldType>
size_t CapacityConstrainedOptimizer<FieldType>::minimize_lagrangian() {
    if (_parameters.inner_solver == "newton") {
        return minimize_lagrangian_by_newton();
    }

    return minimize_lagrangian_by_lbfgs();
}

// Trust-region Newton on the augmented Lagrangian, sharing the step solver
// with the unconstrained relaxation. The curvature is exact for the energy
// and Gauss-Newton for the penalty, which is the term the growing penalty
// makes ill-conditioned for a history-based method.
template<fields::spherical::Field FieldType>
size_t CapacityConstrainedOptimizer<FieldType>::minimize_lagrangian_by_newton() {
    TrustRegionStep solver(
        _parameters.newton.max_conjugate_gradient_iterations,
        _parameters.newton.conjugate_gradient_tolerance
    );

    std::vector<Vector3> current = site_points();
    SphereState state = _lagrangian.sphere_state(*_sphere);
    LagrangianEvaluation evaluation = _lagrangian.evaluate(*_sphere, state, _multipliers, _penalty);
    double radius = _parameters.newton.initial_trust_radius;
    size_t iterations = 0;

    // A rejected step leaves the iterate where it was, so the gradient and
    // the curvature there are still the ones just computed.
    std::vector<Vector3> gradient;
    std::optional<CurvatureOperator> hessian;

    while (iterations < _parameters.max_inner_iterations) {
        if (!hessian.has_value()) {
            gradient = tangential_gradient(evaluation.site_gradients, current);
            hessian.emplace(curvature_operator(current, state, evaluation, gradient));
        }

        if (gradient_norm(gradient) <= _parameters.newton.gradient_tolerance) {
            break;
        }

        TrustRegionStep::Result step = solver.solve(gradient, *hessian, radius);
        ++iterations;

        if (step.predicted_decrease <= 0.0) {
            break;
        }

        std::vector<Vector3> trial = stepped(current, step.step);
        apply_points(trial);
        SphereState trial_state = _lagrangian.sphere_state(*_sphere);
        LagrangianEvaluation trial_evaluation = _lagrangian.evaluate(*_sphere, trial_state, _multipliers, _penalty);
        double ratio = (evaluation.value - trial_evaluation.value) / step.predicted_decrease;

        if (ratio < 0.25) {
            radius *= 0.25;
        } else if (ratio > 0.75 && step.hit_boundary) {
            radius = std::min(2.0 * radius, _parameters.newton.max_trust_radius);
        }

        if (ratio <= _parameters.newton.acceptance_threshold) {
            apply_points(current);

            if (radius < _parameters.newton.minimum_trust_radius) {
                break;
            }

            continue;
        }

        current = std::move(trial);
        state = std::move(trial_state);
        evaluation = std::move(trial_evaluation);
        hessian.reset();
    }

    return iterations;
}

template<fields::spherical::Field FieldType>
typename CapacityConstrainedOptimizer<FieldType>::CurvatureOperator CapacityConstrainedOptimizer<FieldType>::curvature_operator(
    const std::vector<Vector3>& points,
    const SphereState& state,
    const LagrangianEvaluation& evaluation,
    const std::vector<Vector3>& gradient
) const {
    if (_parameters.newton.curvature == "finite-difference") {
        auto gradient_at = [this](const std::vector<Vector3>& displaced) { return tangential_gradient_at(displaced); };
        FiniteDifferenceHessian<decltype(gradient_at)> hessian(points, gradient, gradient_at, FINITE_DIFFERENCE_DISPLACEMENT);
        return {[hessian](const std::vector<Vector3>& directions) { return hessian.multiply(directions); }};
    }

    CapacityConstrainedHessian hessian(
        _curvature.assemble(*_sphere).through_normalization(points, evaluation.site_gradients),
        CapacityJacobian(state, points),
        points,
        _penalty
    );

    return {[hessian](const std::vector<Vector3>& directions) { return hessian.multiply(directions); }};
}

template<fields::spherical::Field FieldType>
std::vector<Vector3> CapacityConstrainedOptimizer<FieldType>::tangential_gradient_at(const std::vector<Vector3>& points) const {
    Sphere sphere;

    for (const Vector3& point : points) {
        sphere.insert(cgal::to_point(VectorS2(point)));
    }

    SphereState state = _lagrangian.sphere_state(sphere);
    LagrangianEvaluation evaluation = _lagrangian.evaluate(sphere, state, _multipliers, _penalty);

    return tangential_gradient(evaluation.site_gradients, points);
}

template<fields::spherical::Field FieldType>
std::vector<Vector3> CapacityConstrainedOptimizer<FieldType>::site_points() const {
    std::vector<Vector3> points;
    points.reserve(_sphere->size());

    for (size_t k = 0; k < _sphere->size(); ++k) {
        points.push_back(to_vector3(_sphere->site(k)));
    }

    return points;
}

template<fields::spherical::Field FieldType>
void CapacityConstrainedOptimizer<FieldType>::apply_points(const std::vector<Vector3>& points) {
    auto sphere = std::make_unique<Sphere>();

    for (const Vector3& point : points) {
        sphere->insert(cgal::to_point(VectorS2(point.normalized())));
    }

    _sphere = std::move(sphere);
}

template<fields::spherical::Field FieldType>
std::vector<Vector3> CapacityConstrainedOptimizer<FieldType>::tangential_gradient(
    const std::vector<Vector3>& site_gradients,
    const std::vector<Vector3>& points
) {
    std::vector<Vector3> gradient(points.size());

    for (size_t k = 0; k < points.size(); ++k) {
        gradient[k] = Normalization(points[k]).gradient(site_gradients[k]);
    }

    return gradient;
}

template<fields::spherical::Field FieldType>
double CapacityConstrainedOptimizer<FieldType>::gradient_norm(const std::vector<Vector3>& gradient) {
    double sum = 0.0;

    for (const Vector3& entry : gradient) {
        sum += entry.squaredNorm();
    }

    return std::sqrt(sum);
}

template<fields::spherical::Field FieldType>
std::vector<Vector3> CapacityConstrainedOptimizer<FieldType>::stepped(
    const std::vector<Vector3>& points,
    const std::vector<Vector3>& step
) {
    std::vector<Vector3> result(points.size());

    for (size_t k = 0; k < points.size(); ++k) {
        result[k] = (points[k] + step[k]).normalized();
    }

    return result;
}

template<fields::spherical::Field FieldType>
size_t CapacityConstrainedOptimizer<FieldType>::minimize_lagrangian_by_lbfgs() {
    LBFGSpp::LBFGSParam<double> parameters;
    parameters.m = static_cast<int>(_parameters.lbfgs_history);
    parameters.epsilon = 1e-12;
    parameters.epsilon_rel = 1e-10;
    parameters.past = 5;
    parameters.delta = 1e-12;
    parameters.max_iterations = static_cast<int>(_parameters.max_inner_iterations);

    LBFGSpp::LBFGSSolver<double> solver(parameters);
    Objective objective(*this);
    Eigen::VectorXd x = sites_to_vector();
    double value = 0.0;
    size_t iterations = 0;

    try {
        iterations = static_cast<size_t>(solver.minimize(objective, x, value));
    } catch (const std::runtime_error&) {
        iterations = objective.evaluations();
    } catch (const std::logic_error&) {
        iterations = objective.evaluations();
    }

    apply_sites(objective.best_x());
    return iterations;
}

template<fields::spherical::Field FieldType>
LagrangianEvaluation CapacityConstrainedOptimizer<FieldType>::evaluate() const {
    return _lagrangian.evaluate(*_sphere, _multipliers, _penalty);
}

template<fields::spherical::Field FieldType>
void CapacityConstrainedOptimizer<FieldType>::update_multipliers(
    const LagrangianEvaluation& evaluation,
    double previous_violation
) {
    for (size_t i = 0; i < _multipliers.size(); ++i) {
        _multipliers[i] += _penalty * evaluation.capacity_errors[i];
    }

    bool insufficient_decrease =
        evaluation.max_absolute_capacity_error() > _parameters.required_violation_decrease * previous_violation;
    bool below_cap = _penalty * _parameters.penalty_growth <= _initial_penalty * _parameters.max_penalty_growth_factor;

    if (insufficient_decrease && below_cap) {
        _penalty *= _parameters.penalty_growth;
    }
}

template<fields::spherical::Field FieldType>
void CapacityConstrainedOptimizer<FieldType>::print_progress(
    size_t outer_iteration,
    size_t inner_iterations,
    const LagrangianEvaluation& evaluation
) const {
    std::cout << "  " << std::setw(8) << std::left << "CCVT" << std::right <<
        std::setw(3) << outer_iteration << " (" << std::setw(4) << inner_iterations << " " << _parameters.inner_solver << ")" <<
        ": CVT energy " << std::scientific << std::setprecision(6) << evaluation.cvt_energy <<
        ", capacity RMS " << std::setprecision(3) << _report.relative_rms_capacity_error <<
        ", penalty " << std::setprecision(2) << _penalty <<
        std::defaultfloat << std::endl;
}

template<fields::spherical::Field FieldType>
Eigen::VectorXd CapacityConstrainedOptimizer<FieldType>::sites_to_vector() const {
    Eigen::VectorXd x(3 * _sphere->size());

    for (size_t i = 0; i < _sphere->size(); ++i) {
        x.segment<3>(3 * static_cast<Eigen::Index>(i)) = to_vector3(_sphere->site(i));
    }

    return x;
}

template<fields::spherical::Field FieldType>
void CapacityConstrainedOptimizer<FieldType>::apply_sites(const Eigen::VectorXd& x) {
    for (size_t i = 0; i < _sphere->size(); ++i) {
        Vector3 site = x.segment<3>(3 * static_cast<Eigen::Index>(i)).normalized();
        _sphere->update_site(i, cgal::to_point(site));
    }
}

template<fields::spherical::Field FieldType>
Eigen::VectorXd CapacityConstrainedOptimizer<FieldType>::chain_rule_gradient(
    const Eigen::VectorXd& x,
    const std::vector<Vector3>& site_gradients
) const {
    Eigen::VectorXd gradient(x.size());

    for (size_t i = 0; i < site_gradients.size(); ++i) {
        Eigen::Index offset = 3 * static_cast<Eigen::Index>(i);
        Normalization normalization(x.segment<3>(offset));
        gradient.segment<3>(offset) = normalization.gradient(site_gradients[i]);
    }

    return gradient;
}

template<fields::spherical::Field FieldType>
CapacityConstrainedOptimizer<FieldType>::Objective::Objective(CapacityConstrainedOptimizer& optimizer) :
    _optimizer(optimizer),
    _best_x(optimizer.sites_to_vector()) {
}

template<fields::spherical::Field FieldType>
double CapacityConstrainedOptimizer<FieldType>::Objective::operator()(
    const Eigen::VectorXd& x,
    Eigen::VectorXd& gradient
) {
    _optimizer.apply_sites(x);
    LagrangianEvaluation evaluation = _optimizer.evaluate();
    gradient = _optimizer.chain_rule_gradient(x, evaluation.site_gradients);
    ++_evaluations;
    _optimizer._callback(*_optimizer._sphere);

    if (evaluation.value < _best_value) {
        _best_value = evaluation.value;
        _best_x = x;
    }

    return evaluation.value;
}

} // namespace globe::voronoi::spherical

#endif //GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_CAPACITY_CONSTRAINED_OPTIMIZER_HPP_
