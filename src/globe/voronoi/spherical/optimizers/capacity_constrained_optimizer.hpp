#ifndef GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_CAPACITY_CONSTRAINED_OPTIMIZER_HPP_
#define GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_CAPACITY_CONSTRAINED_OPTIMIZER_HPP_

#include "../../augmented_lagrangian_loop.hpp"
#include "../../newton_minimizer.hpp"
#include "../../optimizer_parameters.hpp"
#include "capacity_constrained_lagrangian.hpp"
#include "capacity_constrained_hessian.hpp"
#include "capacity_hessian.hpp"
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

    // The ConstrainedProblem interface the shared outer loop drives.
    [[nodiscard]] size_t minimize(const std::vector<double>& multipliers, double penalty);
    [[nodiscard]] LagrangianEvaluation evaluate(
        const std::vector<double>& multipliers,
        double penalty
    ) const;
    [[nodiscard]] double target_mass() const { return _lagrangian.target_mass(); }
    [[nodiscard]] size_t site_count() const { return _sphere->size(); }

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
    CapacityHessian<FieldType> _constraint_curvature;
    CapacityConstrainedParameters _parameters;
    Callback _callback;
    std::vector<double> _multipliers;
    double _penalty = 0.0;
    double _initial_penalty = 0.0;
    CapacityConstrainedReport _report;

    static constexpr double FINITE_DIFFERENCE_DISPLACEMENT = 1e-5;

    size_t minimize_lagrangian();
    [[nodiscard]] LagrangianEvaluation evaluate() const;

    [[nodiscard]] size_t minimize_lagrangian_by_lbfgs();
    [[nodiscard]] size_t minimize_lagrangian_by_newton();

    // The three curvature models share no type, so both directions the
    // descent asks for are carried as functions. A model with no block
    // diagonal to offer preconditions with the identity.
    struct CurvatureOperator {
        using Transform = std::function<std::vector<Vector3>(const std::vector<Vector3>&)>;

        Transform apply;
        Transform inverse_diagonal;

        [[nodiscard]] std::vector<Vector3> multiply(const std::vector<Vector3>& directions) const {
            return apply(directions);
        }

        [[nodiscard]] std::vector<Vector3> precondition(const std::vector<Vector3>& residuals) const {
            return inverse_diagonal(residuals);
        }
    };

    class NewtonLagrangianModel {
     public:
        struct Trial {
            std::unique_ptr<Sphere> sphere;
            std::vector<Vector3> points;
            DiagramState state;
            LagrangianEvaluation evaluation;
        };

        explicit NewtonLagrangianModel(CapacityConstrainedOptimizer& optimizer);

        [[nodiscard]] double value() const { return _evaluation.value; }
        [[nodiscard]] const std::vector<Vector3>& gradient() const { return _gradient; }
        void refresh_curvature();
        [[nodiscard]] const CurvatureOperator& curvature() const { return *_curvature; }
        [[nodiscard]] Trial trial(const std::vector<Vector3>& step) const;
        [[nodiscard]] double trial_value(const Trial& trial) const { return trial.evaluation.value; }
        void accept(Trial&& trial);

     private:
        CapacityConstrainedOptimizer& _optimizer;
        std::vector<Vector3> _points;
        DiagramState _state;
        LagrangianEvaluation _evaluation;
        std::vector<Vector3> _gradient;
        std::optional<CurvatureOperator> _curvature;
    };

    [[nodiscard]] CurvatureOperator curvature_operator(
        const std::vector<Vector3>& points,
        const DiagramState& state,
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
    _constraint_curvature(field),
    _parameters(parameters),
    _callback(std::move(callback)),
    _multipliers(_sphere->size(), 0.0) {
}

template<fields::spherical::Field FieldType>
std::unique_ptr<Sphere> CapacityConstrainedOptimizer<FieldType>::optimize() {
    _report = augmented_lagrangian_loop(*this, _parameters, _parameters.inner_solver);
    return std::move(_sphere);
}

template<fields::spherical::Field FieldType>
size_t CapacityConstrainedOptimizer<FieldType>::minimize(
    const std::vector<double>& multipliers,
    double penalty
) {
    _multipliers = multipliers;
    _penalty = penalty;
    return minimize_lagrangian();
}

template<fields::spherical::Field FieldType>
LagrangianEvaluation CapacityConstrainedOptimizer<FieldType>::evaluate(
    const std::vector<double>& multipliers,
    double penalty
) const {
    return _lagrangian.evaluate(*_sphere, multipliers, penalty);
}

template<fields::spherical::Field FieldType>
size_t CapacityConstrainedOptimizer<FieldType>::minimize_lagrangian() {
    if (_parameters.inner_solver == "newton") {
        return minimize_lagrangian_by_newton();
    }

    return minimize_lagrangian_by_lbfgs();
}

// Trust-region Newton on the augmented Lagrangian through the shared
// descent loop; this model supplies the sphere -- tangential gradients,
// curvature through the normalization, and renormalising rebuilds. The
// curvature is exact for the energy and, by default, for the constraints;
// Gauss-Newton for the penalty.
template<fields::spherical::Field FieldType>
size_t CapacityConstrainedOptimizer<FieldType>::minimize_lagrangian_by_newton() {
    NewtonLagrangianModel model(*this);
    return newton_minimize(model, _parameters.newton, _parameters.max_inner_iterations);
}

template<fields::spherical::Field FieldType>
CapacityConstrainedOptimizer<FieldType>::NewtonLagrangianModel::NewtonLagrangianModel(
    CapacityConstrainedOptimizer& optimizer
) :
    _optimizer(optimizer),
    _points(optimizer.site_points()),
    _state(optimizer._lagrangian.sphere_state(*optimizer._sphere)),
    _evaluation(optimizer._lagrangian.evaluate(*optimizer._sphere, _state, optimizer._multipliers, optimizer._penalty)),
    _gradient(tangential_gradient(_evaluation.site_gradients, _points)) {
}

template<fields::spherical::Field FieldType>
void CapacityConstrainedOptimizer<FieldType>::NewtonLagrangianModel::refresh_curvature() {
    _curvature.emplace(_optimizer.curvature_operator(_points, _state, _evaluation, _gradient));
}

template<fields::spherical::Field FieldType>
typename CapacityConstrainedOptimizer<FieldType>::NewtonLagrangianModel::Trial
CapacityConstrainedOptimizer<FieldType>::NewtonLagrangianModel::trial(const std::vector<Vector3>& step) const {
    Trial trial;
    trial.points = stepped(_points, step);
    trial.sphere = std::make_unique<Sphere>();

    for (const Vector3& point : trial.points) {
        trial.sphere->insert(cgal::to_point(VectorS2(point)));
    }

    trial.state = _optimizer._lagrangian.sphere_state(*trial.sphere);
    trial.evaluation = _optimizer._lagrangian.evaluate(
        *trial.sphere,
        trial.state,
        _optimizer._multipliers,
        _optimizer._penalty
    );

    return trial;
}

template<fields::spherical::Field FieldType>
void CapacityConstrainedOptimizer<FieldType>::NewtonLagrangianModel::accept(Trial&& trial) {
    _optimizer._sphere = std::move(trial.sphere);
    _points = std::move(trial.points);
    _state = std::move(trial.state);
    _evaluation = std::move(trial.evaluation);
    _gradient = tangential_gradient(_evaluation.site_gradients, _points);
    _curvature.reset();
    _optimizer._callback(*_optimizer._sphere);
}

template<fields::spherical::Field FieldType>
typename CapacityConstrainedOptimizer<FieldType>::CurvatureOperator CapacityConstrainedOptimizer<FieldType>::curvature_operator(
    const std::vector<Vector3>& points,
    const DiagramState& state,
    const LagrangianEvaluation& evaluation,
    const std::vector<Vector3>& gradient
) const {
    if (_parameters.newton.curvature == "finite-difference") {
        auto gradient_at = [this](const std::vector<Vector3>& displaced) { return tangential_gradient_at(displaced); };
        FiniteDifferenceHessian<decltype(gradient_at)> hessian(points, gradient, gradient_at, FINITE_DIFFERENCE_DISPLACEMENT);

        return {
            [hessian](const std::vector<Vector3>& directions) { return hessian.multiply(directions); },
            [](const std::vector<Vector3>& residuals) { return residuals; }
        };
    }

    HessianBlocks blocks = _curvature.assemble(*_sphere);

    if (_parameters.newton.curvature == "exact") {
        std::vector<double> weights(_sphere->size());

        for (size_t k = 0; k < weights.size(); ++k) {
            weights[k] = _multipliers[k] + _penalty * evaluation.capacity_errors[k];
        }

        blocks = blocks.plus(_constraint_curvature.assemble(*_sphere, state, points, weights));
    }

    CapacityConstrainedHessian hessian(
        blocks.template through_manifold<Normalization>(points, evaluation.site_gradients),
        CapacityJacobian(state),
        points,
        _penalty
    );

    return {
        [hessian](const std::vector<Vector3>& directions) { return hessian.multiply(directions); },
        [hessian](const std::vector<Vector3>& residuals) { return hessian.precondition(residuals); }
    };
}

template<fields::spherical::Field FieldType>
std::vector<Vector3> CapacityConstrainedOptimizer<FieldType>::tangential_gradient_at(const std::vector<Vector3>& points) const {
    Sphere sphere;

    for (const Vector3& point : points) {
        sphere.insert(cgal::to_point(VectorS2(point)));
    }

    DiagramState state = _lagrangian.sphere_state(sphere);
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
