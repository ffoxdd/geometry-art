#ifndef GEOMETRY_ART_VORONOI_TRUST_REGION_STEP_HPP_
#define GEOMETRY_ART_VORONOI_TRUST_REGION_STEP_HPP_

#include "../types.hpp"
#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <vector>

namespace geometry_art::voronoi {

// A curvature operator, plus the block-diagonal inverse the conjugate
// gradient preconditions with. An operator with no useful block diagonal
// returns its argument, which is the unpreconditioned method.
template<typename T>
concept LinearOperator = requires(const T& operator_, const std::vector<Vector3>& directions) {
    { operator_.multiply(directions) } -> std::convertible_to<std::vector<Vector3>>;
    { operator_.precondition(directions) } -> std::convertible_to<std::vector<Vector3>>;
};

// Steihaug's truncated conjugate gradient: minimises the quadratic model
// g·p + p·Hp/2 over the ball of the given radius, stopping early at the
// boundary or at a direction of negative curvature. Only Hessian-vector
// products are needed, so the Hessian is never factorised.
//
// The iteration is preconditioned by the operator's block diagonal while the
// ball stays Euclidean, so the trust radius keeps its plain meaning as how
// far a site may move.
class TrustRegionStep {
 public:
    struct Result {
        std::vector<Vector3> step;
        double predicted_decrease = 0.0;
        size_t iterations = 0;
        bool hit_boundary = false;
    };

    TrustRegionStep(size_t max_iterations, double relative_tolerance);

    template<LinearOperator OperatorType>
    [[nodiscard]] Result solve(
        const std::vector<Vector3>& gradient,
        const OperatorType& hessian,
        double radius
    ) const;

 private:
    size_t _max_iterations;
    double _relative_tolerance;

    [[nodiscard]] static double dot(const std::vector<Vector3>& left, const std::vector<Vector3>& right);
    [[nodiscard]] static double norm(const std::vector<Vector3>& value);
    [[nodiscard]] static std::vector<Vector3> combine(
        const std::vector<Vector3>& base,
        double scale,
        const std::vector<Vector3>& direction
    );
    [[nodiscard]] static double boundary_scale(
        const std::vector<Vector3>& base,
        const std::vector<Vector3>& direction,
        double radius
    );
};

inline TrustRegionStep::TrustRegionStep(size_t max_iterations, double relative_tolerance) :
    _max_iterations(max_iterations),
    _relative_tolerance(relative_tolerance) {
}

template<LinearOperator OperatorType>
TrustRegionStep::Result TrustRegionStep::solve(
    const std::vector<Vector3>& gradient,
    const OperatorType& hessian,
    double radius
) const {
    Result result;
    result.step.assign(gradient.size(), Vector3::Zero());

    double gradient_norm = norm(gradient);
    double tolerance = std::min(_relative_tolerance, std::sqrt(gradient_norm)) * gradient_norm;

    if (gradient_norm <= tolerance) {
        return result;
    }

    std::vector<Vector3> residual = gradient;
    std::vector<Vector3> preconditioned = hessian.precondition(residual);
    std::vector<Vector3> direction = combine(result.step, -1.0, preconditioned);
    double residual_product = dot(residual, preconditioned);

    while (result.iterations < _max_iterations) {
        std::vector<Vector3> mapped = hessian.multiply(direction);
        double curvature = dot(direction, mapped);
        ++result.iterations;

        if (curvature <= 0.0) {
            result.step = combine(result.step, boundary_scale(result.step, direction, radius), direction);
            result.hit_boundary = true;
            break;
        }

        double step_scale = residual_product / curvature;
        std::vector<Vector3> candidate = combine(result.step, step_scale, direction);

        if (norm(candidate) >= radius) {
            result.step = combine(result.step, boundary_scale(result.step, direction, radius), direction);
            result.hit_boundary = true;
            break;
        }

        result.step = candidate;
        residual = combine(residual, step_scale, mapped);

        if (norm(residual) <= tolerance) {
            break;
        }

        preconditioned = hessian.precondition(residual);
        double next_residual_product = dot(residual, preconditioned);
        double conjugacy = next_residual_product / residual_product;

        for (size_t i = 0; i < direction.size(); ++i) {
            direction[i] = conjugacy * direction[i] - preconditioned[i];
        }

        residual_product = next_residual_product;
    }

    std::vector<Vector3> mapped_step = hessian.multiply(result.step);
    result.predicted_decrease = -(dot(gradient, result.step) + 0.5 * dot(result.step, mapped_step));

    return result;
}

inline double TrustRegionStep::dot(const std::vector<Vector3>& left, const std::vector<Vector3>& right) {
    double sum = 0.0;

    for (size_t i = 0; i < left.size(); ++i) {
        sum += left[i].dot(right[i]);
    }

    return sum;
}

inline double TrustRegionStep::norm(const std::vector<Vector3>& value) {
    return std::sqrt(dot(value, value));
}

inline std::vector<Vector3> TrustRegionStep::combine(
    const std::vector<Vector3>& base,
    double scale,
    const std::vector<Vector3>& direction
) {
    std::vector<Vector3> result(base.size());

    for (size_t i = 0; i < base.size(); ++i) {
        result[i] = base[i] + scale * direction[i];
    }

    return result;
}

// The positive root of ||base + scale * direction|| = radius.
inline double TrustRegionStep::boundary_scale(
    const std::vector<Vector3>& base,
    const std::vector<Vector3>& direction,
    double radius
) {
    double direction_square = dot(direction, direction);

    if (direction_square <= 0.0) {
        return 0.0;
    }

    double cross = dot(base, direction);
    double slack = radius * radius - dot(base, base);
    double discriminant = std::max(0.0, cross * cross + direction_square * slack);

    return (-cross + std::sqrt(discriminant)) / direction_square;
}

} // namespace geometry_art::voronoi

#endif //GEOMETRY_ART_VORONOI_TRUST_REGION_STEP_HPP_
