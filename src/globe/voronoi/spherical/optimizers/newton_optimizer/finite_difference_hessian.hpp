#ifndef GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_NEWTON_OPTIMIZER_FINITE_DIFFERENCE_HESSIAN_HPP_
#define GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_NEWTON_OPTIMIZER_FINITE_DIFFERENCE_HESSIAN_HPP_

#include "../../../../types.hpp"
#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <utility>
#include <vector>

namespace globe::voronoi::spherical {

template<typename T>
concept TangentialGradientFunction = requires(const T& function, const std::vector<Vector3>& points) {
    { function(points) } -> std::convertible_to<std::vector<Vector3>>;
};

// The curvature of a function of site directions, read off a forward
// difference of its gradient along the step path the optimiser uses. Exact
// up to the difference error, including every term an analytic model may
// leave out, which is what makes it the reference the models are held to.
template<TangentialGradientFunction GradientFunction>
class FiniteDifferenceHessian {
 public:
    FiniteDifferenceHessian(
        std::vector<Vector3> points,
        std::vector<Vector3> gradient,
        GradientFunction gradient_at,
        double displacement
    );

    [[nodiscard]] std::vector<Vector3> multiply(const std::vector<Vector3>& directions) const;

 private:
    std::vector<Vector3> _points;
    std::vector<Vector3> _gradient;
    GradientFunction _gradient_at;
    double _displacement;

    [[nodiscard]] std::vector<Vector3> tangential(const std::vector<Vector3>& value) const;
    [[nodiscard]] static double largest_norm(const std::vector<Vector3>& value);
};

template<TangentialGradientFunction GradientFunction>
FiniteDifferenceHessian<GradientFunction>::FiniteDifferenceHessian(
    std::vector<Vector3> points,
    std::vector<Vector3> gradient,
    GradientFunction gradient_at,
    double displacement
) :
    _points(std::move(points)),
    _gradient(std::move(gradient)),
    _gradient_at(std::move(gradient_at)),
    _displacement(displacement) {
}

template<TangentialGradientFunction GradientFunction>
std::vector<Vector3> FiniteDifferenceHessian<GradientFunction>::multiply(const std::vector<Vector3>& directions) const {
    std::vector<Vector3> tangent = tangential(directions);
    double scale = largest_norm(tangent);

    if (scale == 0.0) {
        return std::vector<Vector3>(directions.size(), Vector3::Zero());
    }

    double epsilon = _displacement / scale;
    std::vector<Vector3> displaced(_points.size());

    for (size_t k = 0; k < _points.size(); ++k) {
        displaced[k] = (_points[k] + epsilon * tangent[k]).normalized();
    }

    std::vector<Vector3> displaced_gradient = _gradient_at(displaced);
    std::vector<Vector3> difference(_points.size());

    for (size_t k = 0; k < _points.size(); ++k) {
        difference[k] = (displaced_gradient[k] - _gradient[k]) / epsilon;
    }

    return tangential(difference);
}

template<TangentialGradientFunction GradientFunction>
std::vector<Vector3> FiniteDifferenceHessian<GradientFunction>::tangential(const std::vector<Vector3>& value) const {
    std::vector<Vector3> result(value.size());

    for (size_t k = 0; k < value.size(); ++k) {
        result[k] = value[k] - value[k].dot(_points[k]) * _points[k];
    }

    return result;
}

template<TangentialGradientFunction GradientFunction>
double FiniteDifferenceHessian<GradientFunction>::largest_norm(const std::vector<Vector3>& value) {
    double largest = 0.0;

    for (const Vector3& entry : value) {
        largest = std::max(largest, entry.norm());
    }

    return largest;
}

} // namespace globe::voronoi::spherical

#endif //GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_NEWTON_OPTIMIZER_FINITE_DIFFERENCE_HESSIAN_HPP_
