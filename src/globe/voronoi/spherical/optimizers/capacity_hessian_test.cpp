#include "capacity_hessian.hpp"
#include "capacity_constrained_lagrangian.hpp"
#include "../../capacity_jacobian.hpp"
#include "../../../fields/scalar/noise_field.hpp"
#include "../../../fields/spherical/piecewise_polynomial_field.hpp"
#include "../../../fields/spherical/polynomial_field.hpp"
#include "../../../fields/spherical/powell_sabin_projection.hpp"
#include "../../../generators/spherical/fibonacci_point_generator.hpp"
#include "../../../geometry/spherical/triangle_mesh.hpp"
#include "../../../math/interval.hpp"
#include "../../../testing/macros.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <memory>
#include <vector>

using namespace globe;
using fields::scalar::NoiseField;
using fields::spherical::PiecewisePolynomialField;
using fields::spherical::PolynomialField;
using fields::spherical::PowellSabinProjection;
using geometry::spherical::TriangleMesh;
using voronoi::spherical::CapacityConstrainedLagrangian;
using voronoi::spherical::CapacityHessian;
using voronoi::CapacityJacobian;
using voronoi::HessianBlocks;
using voronoi::spherical::Sphere;
using voronoi::DiagramState;

namespace {

constexpr double DISPLACEMENT = 1e-6;

std::vector<Vector3> sites_on_sphere(size_t count) {
    std::vector<Vector3> sites;

    for (const VectorS2& point : generators::spherical::FibonacciPointGenerator().generate(count)) {
        sites.push_back(point);
    }

    return sites;
}

std::vector<double> arbitrary_weights(size_t count) {
    std::vector<double> weights(count);

    for (size_t k = 0; k < count; ++k) {
        weights[k] = std::sin(3.0 * static_cast<double>(k) + 0.7);
    }

    return weights;
}

std::unique_ptr<Sphere> build_sphere(const std::vector<Vector3>& sites) {
    auto sphere = std::make_unique<Sphere>();

    for (const Vector3& site : sites) {
        sphere->insert(cgal::to_point(VectorS2(site.normalized())));
    }

    return sphere;
}

std::vector<Vector3> tangential(const std::vector<Vector3>& value, const std::vector<Vector3>& points) {
    std::vector<Vector3> result(value.size());

    for (size_t k = 0; k < value.size(); ++k) {
        Vector3 site = points[k].normalized();
        result[k] = value[k] - value[k].dot(site) * site;
    }

    return result;
}

template<fields::spherical::Field FieldType>
std::vector<Vector3> weighted_capacity_gradient(
    const FieldType& field,
    const std::vector<Vector3>& points,
    const std::vector<double>& weights
) {
    auto sphere = build_sphere(points);
    CapacityConstrainedLagrangian<FieldType> lagrangian(field, 0.0);
    DiagramState state = lagrangian.sphere_state(*sphere);

    std::vector<Vector3> normalized(points.size());

    for (size_t k = 0; k < points.size(); ++k) {
        normalized[k] = points[k].normalized();
    }

    return tangential(CapacityJacobian(state).transpose_apply(weights), points);
}

template<fields::spherical::Field FieldType>
void expect_matches_finite_differences(const FieldType& field, size_t site_count, double tolerance) {
    std::vector<Vector3> points = sites_on_sphere(site_count);
    std::vector<double> weights = arbitrary_weights(site_count);

    auto sphere = build_sphere(points);
    CapacityConstrainedLagrangian<FieldType> lagrangian(field, 0.0);
    DiagramState state = lagrangian.sphere_state(*sphere);
    std::vector<Vector3> site_gradients = CapacityJacobian(state).transpose_apply(weights);

    HessianBlocks blocks = CapacityHessian<FieldType>(field)
        .assemble(*sphere, state, points, weights)
        .template through_manifold<Normalization>(points, site_gradients);

    for (int trial = 0; trial < 3; ++trial) {
        std::vector<Vector3> direction(points.size());

        for (size_t k = 0; k < points.size(); ++k) {
            direction[k] = Vector3(
                std::sin(1.3 * static_cast<double>(k) + trial),
                std::cos(2.1 * static_cast<double>(k) - trial),
                std::sin(0.7 * static_cast<double>(k) + 2.0 * trial)
            );
        }

        direction = tangential(direction, points);
        std::vector<Vector3> analytic = blocks.multiply(direction);

        std::vector<Vector3> forward(points.size());
        std::vector<Vector3> backward(points.size());

        for (size_t k = 0; k < points.size(); ++k) {
            forward[k] = (points[k] + DISPLACEMENT * direction[k]).normalized();
            backward[k] = (points[k] - DISPLACEMENT * direction[k]).normalized();
        }

        std::vector<Vector3> ahead = weighted_capacity_gradient(field, forward, weights);
        std::vector<Vector3> behind = weighted_capacity_gradient(field, backward, weights);
        std::vector<Vector3> difference(points.size());

        for (size_t k = 0; k < points.size(); ++k) {
            difference[k] = (ahead[k] - behind[k]) / (2.0 * DISPLACEMENT);
        }

        difference = tangential(difference, points);

        double scale = 0.0;

        for (const Vector3& entry : difference) {
            scale = std::max(scale, entry.norm());
        }

        for (size_t k = 0; k < points.size(); ++k) {
            EXPECT_LT((analytic[k] - difference[k]).norm(), tolerance * scale) <<
                "site " << k << " trial " << trial;
        }
    }
}

} // namespace

// The reference every curvature model is held to: the derivative of the
// exact weighted-capacity gradient along the step path the optimizer uses.
// Matching it means no term is missing -- the rotation of each bisector,
// the change of separation, and the sliding of every Voronoi vertex.
TEST(CapacityHessianTest, MatchesFiniteDifferencesOnAQuadraticField) {
    Eigen::Matrix3d form;
    form << 0.7, 0.2, -0.1,
            0.2, -0.4, 0.3,
           -0.1, 0.3, 0.5;

    PolynomialField field = PolynomialField::quadratic(1.0, Vector3(0.1, -0.2, 0.3), form);
    expect_matches_finite_differences(field, 16, 1e-4);
}

TEST(CapacityHessianTest, MatchesFiniteDifferencesOnAConstantField) {
    expect_matches_finite_differences(PolynomialField::constant(1.0), 12, 1e-4);
}

// The curvature of a true function is symmetric, and every one of the three
// motions must conspire to keep it so: the endpoint terms of one edge match
// the endpoint terms of the edges meeting it.
TEST(CapacityHessianTest, IsSymmetricOnTangentDirections) {
    PolynomialField field = PolynomialField::quadratic(1.0, Vector3(0.1, -0.2, 0.3), 0.2 * Matrix3::Identity());
    std::vector<Vector3> points = sites_on_sphere(14);
    std::vector<double> weights = arbitrary_weights(14);

    auto sphere = build_sphere(points);
    CapacityConstrainedLagrangian<PolynomialField> lagrangian(field, 0.0);
    DiagramState state = lagrangian.sphere_state(*sphere);
    std::vector<Vector3> site_gradients = CapacityJacobian(state).transpose_apply(weights);

    HessianBlocks blocks = CapacityHessian<PolynomialField>(field)
        .assemble(*sphere, state, points, weights)
        .template through_manifold<Normalization>(points, site_gradients);

    std::vector<Vector3> first(points.size());
    std::vector<Vector3> second(points.size());

    for (size_t k = 0; k < points.size(); ++k) {
        first[k] = Vector3(std::sin(1.1 * static_cast<double>(k)), std::cos(0.6 * static_cast<double>(k)), 0.4);
        second[k] = Vector3(0.2, std::sin(2.3 * static_cast<double>(k)), std::cos(1.7 * static_cast<double>(k)));
    }

    first = tangential(first, points);
    second = tangential(second, points);

    std::vector<Vector3> applied_first = blocks.multiply(first);
    std::vector<Vector3> applied_second = blocks.multiply(second);

    double forward = 0.0;
    double backward = 0.0;
    double scale = 0.0;

    for (size_t k = 0; k < points.size(); ++k) {
        forward += second[k].dot(applied_first[k]);
        backward += first[k].dot(applied_second[k]);
        scale = std::max(scale, applied_first[k].norm());
    }

    EXPECT_NEAR(forward, backward, 1e-10 * std::max(scale, 1.0));
}

TEST(CapacityHessianTest, EXPENSIVE_MatchesFiniteDifferencesOnTheSmoothNoiseField) {
    REQUIRE_EXPENSIVE();

    NoiseField noise(Interval(0.2, 1.0));
    PowellSabinProjection projection;
    PiecewisePolynomialField field = projection.project(TriangleMesh::icosphere(2), noise).field;

    expect_matches_finite_differences(field, 20, 1e-3);
}
