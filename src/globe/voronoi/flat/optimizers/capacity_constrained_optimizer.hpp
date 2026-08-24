#ifndef GLOBEART_SRC_GLOBE_VORONOI_FLAT_OPTIMIZERS_CAPACITY_CONSTRAINED_OPTIMIZER_HPP_
#define GLOBEART_SRC_GLOBE_VORONOI_FLAT_OPTIMIZERS_CAPACITY_CONSTRAINED_OPTIMIZER_HPP_

#include "capacity_constrained_lagrangian.hpp"
#include "capacity_hessian.hpp"
#include "cvt_hessian.hpp"
#include "lloyd_optimizer.hpp"
#include "../core/torus.hpp"
#include "../../augmented_lagrangian_loop.hpp"
#include "../../capacity_jacobian.hpp"
#include "../../hessian_blocks.hpp"
#include "../../lagrangian_evaluation.hpp"
#include "../../newton_minimizer.hpp"
#include "../../optimizer_parameters.hpp"
#include "../../state.hpp"
#include "../../../fields/flat/field.hpp"
#include "../../../types.hpp"
#include <CGAL/assertions.h>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

namespace globe::voronoi::flat {

// Capacity-constrained CVT on the flat torus: the shared augmented
// Lagrangian loop around the shared trust-region Newton descent, with this
// class supplying only the geometry -- state assembly, curvature and the
// wrapping rebuild. The curvature model matches the sphere's: exact for the
// energy and, by default, for the constraints; Gauss-Newton for the penalty.
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

    // The ConstrainedProblem interface the shared outer loop drives.
    [[nodiscard]] size_t minimize(const std::vector<double>& multipliers, double penalty);
    [[nodiscard]] LagrangianEvaluation evaluate(const std::vector<double>& multipliers, double penalty) const;
    [[nodiscard]] double target_mass() const { return _lagrangian.target_mass(); }
    [[nodiscard]] size_t site_count() const { return _torus->size(); }

 private:
    struct CurvatureOperator {
        HessianBlocks blocks;
        CapacityJacobian jacobian;
        double penalty;

        [[nodiscard]] std::vector<Vector3> multiply(const std::vector<Vector3>& directions) const;
    };

    // One inner Newton descent at fixed multipliers and penalty.
    class Model {
     public:
        struct Trial {
            std::unique_ptr<Torus> torus;
            DiagramState state;
            LagrangianEvaluation evaluation;
        };

        Model(
            CapacityConstrainedOptimizer& optimizer,
            const std::vector<double>& multipliers,
            double penalty
        );

        [[nodiscard]] double value() const { return _evaluation.value; }
        [[nodiscard]] const std::vector<Vector3>& gradient() const { return _evaluation.site_gradients; }
        void refresh_curvature();
        [[nodiscard]] const CurvatureOperator& curvature() const { return *_curvature; }
        [[nodiscard]] Trial trial(const std::vector<Vector3>& step) const;
        [[nodiscard]] double trial_value(const Trial& trial) const { return trial.evaluation.value; }
        void accept(Trial&& trial);

     private:
        CapacityConstrainedOptimizer& _optimizer;
        const std::vector<double>& _multipliers;
        double _penalty;
        DiagramState _state;
        LagrangianEvaluation _evaluation;
        std::optional<CurvatureOperator> _curvature;
    };

    std::unique_ptr<Torus> _torus;
    CapacityConstrainedLagrangian<FieldType> _lagrangian;
    CvtHessian<FieldType> _curvature;
    CapacityHessian<FieldType> _constraint_curvature;
    CapacityConstrainedParameters _parameters;
    Callback _callback;
    CapacityConstrainedReport _report;

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
    _callback(std::move(callback)) {
    CGAL_precondition(_parameters.inner_solver == "newton");
}

template<fields::flat::Field FieldType>
std::unique_ptr<Torus> CapacityConstrainedOptimizer<FieldType>::optimize() {
    _report = augmented_lagrangian_loop(*this, _parameters, "newton");
    return std::move(_torus);
}

template<fields::flat::Field FieldType>
size_t CapacityConstrainedOptimizer<FieldType>::minimize(
    const std::vector<double>& multipliers,
    double penalty
) {
    Model model(*this, multipliers, penalty);
    return newton_minimize(model, _parameters.newton, _parameters.max_inner_iterations);
}

template<fields::flat::Field FieldType>
LagrangianEvaluation CapacityConstrainedOptimizer<FieldType>::evaluate(
    const std::vector<double>& multipliers,
    double penalty
) const {
    return _lagrangian.evaluate(*_torus, multipliers, penalty);
}

template<fields::flat::Field FieldType>
CapacityConstrainedOptimizer<FieldType>::Model::Model(
    CapacityConstrainedOptimizer& optimizer,
    const std::vector<double>& multipliers,
    double penalty
) :
    _optimizer(optimizer),
    _multipliers(multipliers),
    _penalty(penalty),
    _state(optimizer._lagrangian.diagram_state(*optimizer._torus)),
    _evaluation(optimizer._lagrangian.evaluate(*optimizer._torus, _state, multipliers, penalty)) {
}

template<fields::flat::Field FieldType>
void CapacityConstrainedOptimizer<FieldType>::Model::refresh_curvature() {
    HessianBlocks blocks = _optimizer._curvature.assemble(*_optimizer._torus);

    if (_optimizer._parameters.newton.curvature == "exact") {
        std::vector<double> weights(_optimizer._torus->size());

        for (size_t k = 0; k < weights.size(); ++k) {
            weights[k] = _multipliers[k] + _penalty * _evaluation.capacity_errors[k];
        }

        blocks = blocks.plus(_optimizer._constraint_curvature.assemble(*_optimizer._torus, weights));
    }

    _curvature.emplace(CurvatureOperator{std::move(blocks), CapacityJacobian(_state), _penalty});
}

template<fields::flat::Field FieldType>
typename CapacityConstrainedOptimizer<FieldType>::Model::Trial
CapacityConstrainedOptimizer<FieldType>::Model::trial(const std::vector<Vector3>& step) const {
    std::vector<Vector3> points;
    points.reserve(_optimizer._torus->size());

    for (size_t k = 0; k < _optimizer._torus->size(); ++k) {
        points.push_back(_optimizer._torus->site_vector(k) + step[k]);
    }

    Trial trial;
    trial.torus = _optimizer._torus->rebuilt(points);
    trial.state = _optimizer._lagrangian.diagram_state(*trial.torus);
    trial.evaluation = _optimizer._lagrangian.evaluate(*trial.torus, trial.state, _multipliers, _penalty);
    return trial;
}

template<fields::flat::Field FieldType>
void CapacityConstrainedOptimizer<FieldType>::Model::accept(Trial&& trial) {
    _optimizer._torus = std::move(trial.torus);
    _state = std::move(trial.state);
    _evaluation = std::move(trial.evaluation);
    _curvature.reset();
    _optimizer._callback(*_optimizer._torus);
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

} // namespace globe::voronoi::flat

#endif //GLOBEART_SRC_GLOBE_VORONOI_FLAT_OPTIMIZERS_CAPACITY_CONSTRAINED_OPTIMIZER_HPP_
