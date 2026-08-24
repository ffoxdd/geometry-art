#ifndef GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_NEWTON_OPTIMIZER_NEWTON_OPTIMIZER_HPP_
#define GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_NEWTON_OPTIMIZER_NEWTON_OPTIMIZER_HPP_

#include "trust_region_step.hpp"
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

    std::unique_ptr<Sphere> _sphere;
    CapacityConstrainedLagrangian<FieldType> _lagrangian;
    CvtHessian<FieldType> _hessian;
    NewtonParameters _parameters;
    Callback _callback;
    NewtonReport _report;

    [[nodiscard]] std::vector<Vector3> points() const;
    void apply(const std::vector<Vector3>& points);
    [[nodiscard]] Evaluation evaluate() const;
    [[nodiscard]] double trust_radius_after(double radius, double ratio, bool hit_boundary) const;

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
    TrustRegionStep solver(
        _parameters.max_conjugate_gradient_iterations,
        _parameters.conjugate_gradient_tolerance
    );

    std::vector<Vector3> current = points();
    Evaluation evaluation = evaluate();
    double radius = _parameters.initial_trust_radius;

    // A rejected step leaves the iterate where it was, so the curvature
    // there is still the one just assembled.
    std::optional<HessianBlocks> blocks;

    while (_report.iterations < _parameters.max_iterations) {
        _report.gradient_norm = norm(evaluation.gradient);

        if (_report.gradient_norm <= _parameters.gradient_tolerance) {
            _report.converged = true;
            break;
        }

        if (!blocks.has_value()) {
            blocks = _hessian.assemble(*_sphere).through_normalization(current, evaluation.site_gradients);
        }

        TrustRegionStep::Result step = solver.solve(evaluation.gradient, *blocks, radius);
        ++_report.iterations;

        if (step.predicted_decrease <= 0.0) {
            _report.stalled = true;
            break;
        }

        std::vector<Vector3> trial = advanced(current, step.step);
        apply(trial);
        Evaluation trial_evaluation = evaluate();
        double ratio = (evaluation.energy - trial_evaluation.energy) / step.predicted_decrease;
        radius = trust_radius_after(radius, ratio, step.hit_boundary);

        if (ratio <= _parameters.acceptance_threshold) {
            apply(current);

            if (radius < _parameters.minimum_trust_radius) {
                _report.stalled = true;
                break;
            }

            continue;
        }

        current = std::move(trial);
        evaluation = std::move(trial_evaluation);
        blocks.reset();
        ++_report.accepted_steps;
        _callback(*_sphere);
    }

    _report.cvt_energy = evaluation.energy;
    return std::move(_sphere);
}

template<fields::spherical::Field FieldType>
double NewtonOptimizer<FieldType>::trust_radius_after(double radius, double ratio, bool hit_boundary) const {
    if (ratio < 0.25) {
        return 0.25 * radius;
    }

    if (ratio > 0.75 && hit_boundary) {
        return std::min(2.0 * radius, _parameters.max_trust_radius);
    }

    return radius;
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

// Rebuilt rather than updated site by site: bulk insertion is faster, and
// insertion order is what indexes the cells.
template<fields::spherical::Field FieldType>
void NewtonOptimizer<FieldType>::apply(const std::vector<Vector3>& points) {
    auto sphere = std::make_unique<Sphere>();

    for (const Vector3& point : points) {
        sphere->insert(cgal::to_point(VectorS2(point.normalized())));
    }

    _sphere = std::move(sphere);
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
