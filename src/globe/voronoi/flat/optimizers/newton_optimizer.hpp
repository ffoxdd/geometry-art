#ifndef GLOBEART_SRC_GLOBE_VORONOI_FLAT_OPTIMIZERS_NEWTON_OPTIMIZER_HPP_
#define GLOBEART_SRC_GLOBE_VORONOI_FLAT_OPTIMIZERS_NEWTON_OPTIMIZER_HPP_

#include "capacity_constrained_lagrangian.hpp"
#include "cvt_hessian.hpp"
#include "lloyd_optimizer.hpp"
#include "../core/torus.hpp"
#include "../../hessian_blocks.hpp"
#include "../../lagrangian_evaluation.hpp"
#include "../../newton_minimizer.hpp"
#include "../../optimizer_parameters.hpp"
#include "../../../fields/flat/field.hpp"
#include "../../../types.hpp"
#include <cmath>
#include <cstddef>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

namespace globe::voronoi::flat {

// Unconstrained CVT relaxation on the flat torus through the shared
// trust-region descent: the Lagrangian at zero multipliers and penalty is
// the CVT energy, and its exact curvature is the CvtHessian alone.
template<fields::flat::Field FieldType>
class NewtonOptimizer {
 public:
    NewtonOptimizer(
        std::unique_ptr<Torus> torus,
        FieldType field,
        NewtonParameters parameters,
        Callback callback
    );

    std::unique_ptr<Torus> optimize();
    [[nodiscard]] const NewtonReport& report() const { return _report; }

 private:
    class Model {
     public:
        struct Trial {
            std::unique_ptr<Torus> torus;
            LagrangianEvaluation evaluation;
        };

        explicit Model(NewtonOptimizer& optimizer);

        [[nodiscard]] double value() const { return _evaluation.cvt_energy; }
        [[nodiscard]] const std::vector<Vector3>& gradient() const { return _evaluation.site_gradients; }
        void refresh_curvature();
        [[nodiscard]] const HessianBlocks& curvature() const { return *_curvature; }
        [[nodiscard]] Trial trial(const std::vector<Vector3>& step) const;
        [[nodiscard]] double trial_value(const Trial& trial) const { return trial.evaluation.cvt_energy; }
        void accept(Trial&& trial);

        [[nodiscard]] const LagrangianEvaluation& evaluation() const { return _evaluation; }
        [[nodiscard]] size_t accepted_steps() const { return _accepted_steps; }

     private:
        NewtonOptimizer& _optimizer;
        LagrangianEvaluation _evaluation;
        std::optional<HessianBlocks> _curvature;
        size_t _accepted_steps = 0;
    };

    std::unique_ptr<Torus> _torus;
    CapacityConstrainedLagrangian<FieldType> _lagrangian;
    CvtHessian<FieldType> _hessian;
    NewtonParameters _parameters;
    Callback _callback;
    NewtonReport _report;

    [[nodiscard]] LagrangianEvaluation evaluate() const;
};

template<fields::flat::Field FieldType>
NewtonOptimizer<FieldType>::NewtonOptimizer(
    std::unique_ptr<Torus> torus,
    FieldType field,
    NewtonParameters parameters,
    Callback callback
) :
    _torus(std::move(torus)),
    _lagrangian(field, field.total_mass() / static_cast<double>(_torus->size())),
    _hessian(field),
    _parameters(parameters),
    _callback(std::move(callback)) {
}

template<fields::flat::Field FieldType>
std::unique_ptr<Torus> NewtonOptimizer<FieldType>::optimize() {
    Model model(*this);
    _report.iterations = newton_minimize(model, _parameters, _parameters.max_iterations);
    _report.accepted_steps = model.accepted_steps();
    _report.cvt_energy = model.evaluation().cvt_energy;

    double gradient_norm = 0.0;

    for (const Vector3& entry : model.evaluation().site_gradients) {
        gradient_norm += entry.squaredNorm();
    }

    _report.gradient_norm = std::sqrt(gradient_norm);
    _report.converged = _report.gradient_norm <= _parameters.gradient_tolerance;
    _report.stalled = !_report.converged && _report.iterations < _parameters.max_iterations;

    return std::move(_torus);
}

template<fields::flat::Field FieldType>
NewtonOptimizer<FieldType>::Model::Model(NewtonOptimizer& optimizer) :
    _optimizer(optimizer),
    _evaluation(optimizer.evaluate()) {
}

template<fields::flat::Field FieldType>
LagrangianEvaluation NewtonOptimizer<FieldType>::evaluate() const {
    return _lagrangian.evaluate(*_torus, std::vector<double>(_torus->size(), 0.0), 0.0);
}

template<fields::flat::Field FieldType>
void NewtonOptimizer<FieldType>::Model::refresh_curvature() {
    _curvature.emplace(_optimizer._hessian.assemble(*_optimizer._torus));
}

template<fields::flat::Field FieldType>
typename NewtonOptimizer<FieldType>::Model::Trial
NewtonOptimizer<FieldType>::Model::trial(const std::vector<Vector3>& step) const {
    std::vector<Vector3> points;
    points.reserve(_optimizer._torus->size());

    for (size_t k = 0; k < _optimizer._torus->size(); ++k) {
        points.push_back(_optimizer._torus->site_vector(k) + step[k]);
    }

    Trial trial;
    trial.torus = _optimizer._torus->rebuilt(points);

    std::swap(_optimizer._torus, trial.torus);
    trial.evaluation = _optimizer.evaluate();
    std::swap(_optimizer._torus, trial.torus);

    return trial;
}

template<fields::flat::Field FieldType>
void NewtonOptimizer<FieldType>::Model::accept(Trial&& trial) {
    _optimizer._torus = std::move(trial.torus);
    _evaluation = std::move(trial.evaluation);
    _curvature.reset();
    ++_accepted_steps;
    _optimizer._callback(*_optimizer._torus);
}

} // namespace globe::voronoi::flat

#endif //GLOBEART_SRC_GLOBE_VORONOI_FLAT_OPTIMIZERS_NEWTON_OPTIMIZER_HPP_
