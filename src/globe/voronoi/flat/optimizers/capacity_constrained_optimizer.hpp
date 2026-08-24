#ifndef GLOBEART_SRC_GLOBE_VORONOI_FLAT_OPTIMIZERS_CAPACITY_CONSTRAINED_OPTIMIZER_HPP_
#define GLOBEART_SRC_GLOBE_VORONOI_FLAT_OPTIMIZERS_CAPACITY_CONSTRAINED_OPTIMIZER_HPP_

#include "capacity_constrained_lagrangian.hpp"
#include "capacity_hessian.hpp"
#include "cvt_hessian.hpp"
#include "lloyd_optimizer.hpp"
#include "../core/torus.hpp"
#include "../../capacity_jacobian.hpp"
#include "../../hessian_blocks.hpp"
#include "../../lagrangian_evaluation.hpp"
#include "../../optimizer_parameters.hpp"
#include "../../state.hpp"
#include "../../trust_region_step.hpp"
#include "../../../fields/flat/field.hpp"
#include "../../../types.hpp"
#include <CGAL/assertions.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

namespace globe::voronoi::flat {

// The augmented Lagrangian loop on the flat torus, with a trust-region
// Newton inner solver. The curvature model matches the sphere's: exact for
// the energy, Gauss-Newton for the penalty, and -- by default -- the
// constraints' own second derivatives weighted by multiplier and violation,
// which do not vanish at a solution wherever the density varies.
template<fields::flat::Field FieldType>
class CapacityConstrainedOptimizer {
 public:
    CapacityConstrainedOptimizer(
        std::unique_ptr<Torus> torus,
        FieldType field,
        CapacityConstrainedParameters parameters,
        Callback callback
    );

    std::unique_ptr<Torus> optimize();
    [[nodiscard]] const CapacityConstrainedReport& report() const { return _report; }

 private:
    struct CurvatureOperator {
        HessianBlocks blocks;
        CapacityJacobian jacobian;
        double penalty;

        [[nodiscard]] std::vector<Vector3> multiply(const std::vector<Vector3>& directions) const;
    };

    std::unique_ptr<Torus> _torus;
    CapacityConstrainedLagrangian<FieldType> _lagrangian;
    CvtHessian<FieldType> _curvature;
    CapacityHessian<FieldType> _constraint_curvature;
    CapacityConstrainedParameters _parameters;
    Callback _callback;
    std::vector<double> _multipliers;
    double _penalty = 0.0;
    double _initial_penalty = 0.0;
    CapacityConstrainedReport _report;

    [[nodiscard]] double initial_penalty() const;
    [[nodiscard]] size_t minimize_lagrangian();
    [[nodiscard]] LagrangianEvaluation evaluate() const;
    void update_multipliers(const LagrangianEvaluation& evaluation, double previous_violation);
    void print_progress(size_t outer_iteration, size_t inner_iterations, const LagrangianEvaluation& evaluation) const;

    [[nodiscard]] std::vector<Vector3> site_points() const;
    [[nodiscard]] static double gradient_norm(const std::vector<Vector3>& gradient);
    [[nodiscard]] static std::vector<Vector3> stepped(
        const std::vector<Vector3>& points,
        const std::vector<Vector3>& step
    );
};

template<fields::flat::Field FieldType>
CapacityConstrainedOptimizer<FieldType>::CapacityConstrainedOptimizer(
    std::unique_ptr<Torus> torus,
    FieldType field,
    CapacityConstrainedParameters parameters,
    Callback callback
) :
    _torus(std::move(torus)),
    _lagrangian(field, field.total_mass() / static_cast<double>(_torus->size())),
    _curvature(field),
    _constraint_curvature(field),
    _parameters(parameters),
    _callback(std::move(callback)),
    _multipliers(_torus->size(), 0.0) {
    CGAL_precondition(_parameters.inner_solver == "newton");
}

template<fields::flat::Field FieldType>
std::unique_ptr<Torus> CapacityConstrainedOptimizer<FieldType>::optimize() {
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
        _report.relative_rms_capacity_error =
            evaluation.root_mean_square_capacity_error() / _lagrangian.target_mass();
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

    return std::move(_torus);
}

template<fields::flat::Field FieldType>
double CapacityConstrainedOptimizer<FieldType>::initial_penalty() const {
    LagrangianEvaluation evaluation = evaluate();
    double squared_violation = 0.0;

    for (double error : evaluation.capacity_errors) {
        squared_violation += error * error;
    }

    double target = _lagrangian.target_mass();
    double floor = static_cast<double>(_torus->size()) *
        std::pow(_parameters.relative_capacity_tolerance * target, 2);

    return 2.0 * std::max(evaluation.cvt_energy, std::numeric_limits<double>::min()) /
        std::max(squared_violation, floor);
}

template<fields::flat::Field FieldType>
size_t CapacityConstrainedOptimizer<FieldType>::minimize_lagrangian() {
    TrustRegionStep solver(
        _parameters.newton.max_conjugate_gradient_iterations,
        _parameters.newton.conjugate_gradient_tolerance
    );

    std::vector<Vector3> current = site_points();
    DiagramState state = _lagrangian.diagram_state(*_torus);
    LagrangianEvaluation evaluation = _lagrangian.evaluate(*_torus, state, _multipliers, _penalty);
    double radius = _parameters.newton.initial_trust_radius;
    size_t iterations = 0;

    // A rejected step leaves the iterate where it was, so the gradient and
    // the curvature there are still the ones just computed.
    std::optional<CurvatureOperator> hessian;

    while (iterations < _parameters.max_inner_iterations) {
        if (!hessian.has_value()) {
            HessianBlocks blocks = _curvature.assemble(*_torus);

            if (_parameters.newton.curvature == "exact") {
                std::vector<double> weights(_torus->size());

                for (size_t k = 0; k < weights.size(); ++k) {
                    weights[k] = _multipliers[k] + _penalty * evaluation.capacity_errors[k];
                }

                blocks = blocks.plus(_constraint_curvature.assemble(*_torus, weights));
            }

            hessian.emplace(CurvatureOperator{std::move(blocks), CapacityJacobian(state), _penalty});
        }

        if (gradient_norm(evaluation.site_gradients) <= _parameters.newton.gradient_tolerance) {
            break;
        }

        TrustRegionStep::Result step = solver.solve(evaluation.site_gradients, *hessian, radius);
        ++iterations;

        if (step.predicted_decrease <= 0.0) {
            break;
        }

        std::vector<Vector3> trial_points = stepped(current, step.step);
        auto trial = _torus->rebuilt(trial_points);
        DiagramState trial_state = _lagrangian.diagram_state(*trial);
        LagrangianEvaluation trial_evaluation = _lagrangian.evaluate(*trial, trial_state, _multipliers, _penalty);
        double ratio = (evaluation.value - trial_evaluation.value) / step.predicted_decrease;

        if (ratio < 0.25) {
            radius *= 0.25;
        } else if (ratio > 0.75 && step.hit_boundary) {
            radius = std::min(2.0 * radius, _parameters.newton.max_trust_radius);
        }

        if (ratio <= _parameters.newton.acceptance_threshold) {
            if (radius < _parameters.newton.minimum_trust_radius) {
                break;
            }

            continue;
        }

        _torus = std::move(trial);
        current = site_points();
        state = std::move(trial_state);
        evaluation = std::move(trial_evaluation);
        hessian.reset();
        _callback(*_torus);
    }

    return iterations;
}

template<fields::flat::Field FieldType>
std::vector<Vector3> CapacityConstrainedOptimizer<FieldType>::CurvatureOperator::multiply(
    const std::vector<Vector3>& directions
) const {
    std::vector<Vector3> result = blocks.multiply(directions);

    if (penalty <= 0.0) {
        return result;
    }

    std::vector<Vector3> penalty_term = jacobian.transpose_apply(jacobian.apply(directions));

    for (size_t k = 0; k < result.size(); ++k) {
        result[k] += penalty * penalty_term[k];
    }

    return result;
}

template<fields::flat::Field FieldType>
LagrangianEvaluation CapacityConstrainedOptimizer<FieldType>::evaluate() const {
    return _lagrangian.evaluate(*_torus, _multipliers, _penalty);
}

template<fields::flat::Field FieldType>
void CapacityConstrainedOptimizer<FieldType>::update_multipliers(
    const LagrangianEvaluation& evaluation,
    double previous_violation
) {
    for (size_t i = 0; i < _multipliers.size(); ++i) {
        _multipliers[i] += _penalty * evaluation.capacity_errors[i];
    }

    bool insufficient_decrease =
        evaluation.max_absolute_capacity_error() > _parameters.required_violation_decrease * previous_violation;
    bool below_cap =
        _penalty * _parameters.penalty_growth <= _initial_penalty * _parameters.max_penalty_growth_factor;

    if (insufficient_decrease && below_cap) {
        _penalty *= _parameters.penalty_growth;
    }
}

template<fields::flat::Field FieldType>
void CapacityConstrainedOptimizer<FieldType>::print_progress(
    size_t outer_iteration,
    size_t inner_iterations,
    const LagrangianEvaluation& evaluation
) const {
    std::cout << "  " << std::setw(8) << std::left << "CCVT" << std::right <<
        std::setw(3) << outer_iteration << " (" << std::setw(4) << inner_iterations << " newton)" <<
        ": CVT energy " << std::scientific << std::setprecision(6) << evaluation.cvt_energy <<
        ", capacity RMS " << std::setprecision(3) << _report.relative_rms_capacity_error <<
        ", penalty " << std::setprecision(2) << _penalty <<
        std::defaultfloat << std::endl;
}

template<fields::flat::Field FieldType>
std::vector<Vector3> CapacityConstrainedOptimizer<FieldType>::site_points() const {
    std::vector<Vector3> points;
    points.reserve(_torus->size());

    for (size_t k = 0; k < _torus->size(); ++k) {
        points.push_back(_torus->site_vector(k));
    }

    return points;
}

template<fields::flat::Field FieldType>
double CapacityConstrainedOptimizer<FieldType>::gradient_norm(const std::vector<Vector3>& gradient) {
    double sum = 0.0;

    for (const Vector3& entry : gradient) {
        sum += entry.squaredNorm();
    }

    return std::sqrt(sum);
}

template<fields::flat::Field FieldType>
std::vector<Vector3> CapacityConstrainedOptimizer<FieldType>::stepped(
    const std::vector<Vector3>& points,
    const std::vector<Vector3>& step
) {
    std::vector<Vector3> result(points.size());

    for (size_t k = 0; k < points.size(); ++k) {
        result[k] = points[k] + step[k];
    }

    return result;
}

} // namespace globe::voronoi::flat

#endif //GLOBEART_SRC_GLOBE_VORONOI_FLAT_OPTIMIZERS_CAPACITY_CONSTRAINED_OPTIMIZER_HPP_
