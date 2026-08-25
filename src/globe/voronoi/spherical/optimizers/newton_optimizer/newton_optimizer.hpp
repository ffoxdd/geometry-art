#ifndef GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_NEWTON_OPTIMIZER_NEWTON_OPTIMIZER_HPP_
#define GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_NEWTON_OPTIMIZER_NEWTON_OPTIMIZER_HPP_

#include "../../../hessian_blocks.hpp"
#include "../../../newton_minimizer.hpp"
#include "../../../optimizer_parameters.hpp"
#include "../../../trust_region_step.hpp"
#include "../capacity_constrained_lagrangian.hpp"
#include "../cvt_hessian.hpp"
#include "../../../../types.hpp"
#include "../../../../fields/spherical/field.hpp"
#include "../../../../fields/spherical/polynomial_field.hpp"
#include "../../../../math/normalization.hpp"
#include "../../core/callback.hpp"
#include "../../core/sphere.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace globe::voronoi::spherical {

// Minimises the CVT energy by trust-region Newton over the sites as ambient
// vectors. The energy is scale invariant in each site, so a step is taken in
// the ambient space and the result renormalised.
template<fields::spherical::Field FieldType = fields::spherical::PolynomialField>
class NewtonOptimizer {
 public:
    NewtonOptimizer(
        std::unique_ptr<Sphere> sphere,
        FieldType field,
        NewtonParameters parameters,
        Callback callback
    );

    std::unique_ptr<Sphere> optimize();
    [[nodiscard]] const NewtonReport& report() const { return _report; }

 private:
    struct Evaluation {
        double energy = 0.0;
        std::vector<Vector3> site_gradients;
        std::vector<Vector3> gradient;
    };

    class Model {
     public:
        struct Trial {
            std::unique_ptr<Sphere> sphere;
            std::vector<Vector3> points;
            Evaluation evaluation;
        };

        explicit Model(NewtonOptimizer& optimizer);

        [[nodiscard]] double value() const { return _evaluation.energy; }
        [[nodiscard]] const std::vector<Vector3>& gradient() const { return _evaluation.gradient; }
        void refresh_curvature();
        [[nodiscard]] const PreconditionedBlocks& curvature() const { return *_curvature; }
        [[nodiscard]] Trial trial(const std::vector<Vector3>& step) const;
        [[nodiscard]] double trial_value(const Trial& trial) const { return trial.evaluation.energy; }
        void accept(Trial&& trial);

        [[nodiscard]] const Evaluation& evaluation() const { return _evaluation; }
        [[nodiscard]] size_t accepted_steps() const { return _accepted_steps; }

     private:
        NewtonOptimizer& _optimizer;
        std::vector<Vector3> _points;
        Evaluation _evaluation;
        std::optional<PreconditionedBlocks> _curvature;
        size_t _accepted_steps = 0;
    };

    std::unique_ptr<Sphere> _sphere;
    CapacityConstrainedLagrangian<FieldType> _lagrangian;
    CvtHessian<FieldType> _hessian;
    NewtonParameters _parameters;
    Callback _callback;
    NewtonReport _report;

    [[nodiscard]] std::vector<Vector3> points() const;
    [[nodiscard]] Evaluation evaluate() const;

    [[nodiscard]] static double norm(const std::vector<Vector3>& value);
    [[nodiscard]] static std::vector<Vector3> advanced(
        const std::vector<Vector3>& points,
        const std::vector<Vector3>& step
    );
};

template<fields::spherical::Field FieldType>
NewtonOptimizer<FieldType>::NewtonOptimizer(
    std::unique_ptr<Sphere> sphere,
    FieldType field,
    NewtonParameters parameters,
    Callback callback
) :
    _sphere(std::move(sphere)),
    _lagrangian(field, field.total_mass() / static_cast<double>(_sphere->size())),
    _hessian(field),
    _parameters(parameters),
    _callback(std::move(callback)) {
}

template<fields::spherical::Field FieldType>
std::unique_ptr<Sphere> NewtonOptimizer<FieldType>::optimize() {
    Model model(*this);
    _report.iterations = newton_minimize(model, _parameters, _parameters.max_iterations);
    _report.accepted_steps = model.accepted_steps();
    _report.cvt_energy = model.evaluation().energy;

    double gradient_norm = norm(model.evaluation().gradient);
    _report.gradient_norm = gradient_norm;
    _report.converged = gradient_norm <= _parameters.gradient_tolerance;
    _report.stalled = !_report.converged && _report.iterations < _parameters.max_iterations;

    return std::move(_sphere);
}

template<fields::spherical::Field FieldType>
NewtonOptimizer<FieldType>::Model::Model(NewtonOptimizer& optimizer) :
    _optimizer(optimizer),
    _points(optimizer.points()),
    _evaluation(optimizer.evaluate()) {
}

template<fields::spherical::Field FieldType>
void NewtonOptimizer<FieldType>::Model::refresh_curvature() {
    _curvature.emplace(PreconditionedBlocks::of(
        _optimizer._hessian.assemble(*_optimizer._sphere)
            .template through_manifold<Normalization>(_points, _evaluation.site_gradients)
    ));
}

template<fields::spherical::Field FieldType>
typename NewtonOptimizer<FieldType>::Model::Trial
NewtonOptimizer<FieldType>::Model::trial(const std::vector<Vector3>& step) const {
    Trial trial;
    trial.points = advanced(_points, step);
    trial.sphere = std::make_unique<Sphere>();

    for (const Vector3& point : trial.points) {
        trial.sphere->insert(cgal::to_point(VectorS2(point)));
    }

    std::swap(_optimizer._sphere, trial.sphere);
    trial.evaluation = _optimizer.evaluate();
    std::swap(_optimizer._sphere, trial.sphere);

    return trial;
}

template<fields::spherical::Field FieldType>
void NewtonOptimizer<FieldType>::Model::accept(Trial&& trial) {
    _optimizer._sphere = std::move(trial.sphere);
    _points = std::move(trial.points);
    _evaluation = std::move(trial.evaluation);
    _curvature.reset();
    ++_accepted_steps;
    _optimizer._callback(*_optimizer._sphere);
}

template<fields::spherical::Field FieldType>
std::vector<Vector3> NewtonOptimizer<FieldType>::points() const {
    std::vector<Vector3> result;
    result.reserve(_sphere->size());

    for (size_t k = 0; k < _sphere->size(); ++k) {
        result.push_back(to_vector3(_sphere->site(k)));
    }

    return result;
}

template<fields::spherical::Field FieldType>
typename NewtonOptimizer<FieldType>::Evaluation NewtonOptimizer<FieldType>::evaluate() const {
    LagrangianEvaluation lagrangian = _lagrangian.evaluate(
        *_sphere,
        std::vector<double>(_sphere->size(), 0.0),
        0.0
    );

    Evaluation evaluation;
    evaluation.energy = lagrangian.cvt_energy;
    evaluation.site_gradients = std::move(lagrangian.site_gradients);
    evaluation.gradient.reserve(_sphere->size());

    for (size_t k = 0; k < _sphere->size(); ++k) {
        Normalization normalization(to_vector3(_sphere->site(k)));
        evaluation.gradient.push_back(normalization.gradient(evaluation.site_gradients[k]));
    }

    return evaluation;
}

template<fields::spherical::Field FieldType>
double NewtonOptimizer<FieldType>::norm(const std::vector<Vector3>& value) {
    double sum = 0.0;

    for (const Vector3& entry : value) {
        sum += entry.squaredNorm();
    }

    return std::sqrt(sum);
}

template<fields::spherical::Field FieldType>
std::vector<Vector3> NewtonOptimizer<FieldType>::advanced(
    const std::vector<Vector3>& points,
    const std::vector<Vector3>& step
) {
    std::vector<Vector3> result(points.size());

    for (size_t k = 0; k < points.size(); ++k) {
        result[k] = (points[k] + step[k]).normalized();
    }

    return result;
}

} // namespace globe::voronoi::spherical

#endif //GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_NEWTON_OPTIMIZER_NEWTON_OPTIMIZER_HPP_
