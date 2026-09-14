#include "capacity_hessian.hpp"
#include "capacity_constrained_lagrangian.hpp"
#include "../core/diagram.hpp"
#include "../../capacity_jacobian.hpp"
#include "../../hessian_blocks.hpp"
#include "../../../fields/flat/noise_field.hpp"
#include "../../../fields/flat/piecewise_polynomial_field.hpp"
#include "../../../fields/flat/polynomial_field.hpp"
#include "../../../geometry/planar/domain.hpp"
#include "../../../math/interval.hpp"
#include "../../../testing/flat_scatter.hpp"
#include "../../../testing/macros.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

using namespace geometry_art;
using namespace geometry_art::voronoi;
using namespace geometry_art::voronoi::flat;
using geometry::planar::Domain;
using geometry_art::testing::diagram_of;
using geometry_art::testing::flat_scatter;
using fields::flat::NoiseField;
using fields::flat::PiecewisePolynomialField;
using fields::flat::PolynomialField;
using geometry_art::math::Interval;

namespace {

constexpr double DISPLACEMENT = 1e-6;

std::vector<double> arbitrary_weights(size_t count) {
    std::vector<double> weights(count);

    for (size_t k = 0; k < count; ++k) {
        weights[k] = std::sin(3.0 * static_cast<double>(k) + 0.7);
    }

    return weights;
}

template<fields::flat::Field FieldType>
std::vector<Vector3> weighted_capacity_gradient(
    const FieldType& field,
    const std::vector<Vector2>& sites,
    const Domain& domain,
    const std::vector<double>& weights
) {
    CapacityConstrainedLagrangian<FieldType> lagrangian(field, 0.0);
    auto diagram = diagram_of(sites, domain);

    return CapacityJacobian(lagrangian.diagram_state(*diagram)).transpose_apply(weights);
}

template<fields::flat::Field FieldType>
void expect_matches_finite_differences(
    const FieldType& field,
    size_t count,
    const Domain& domain,
    double tolerance
) {
    std::vector<Vector2> sites = flat_scatter(count, domain.width, domain.height);
    std::vector<double> weights = arbitrary_weights(count);
    auto diagram = diagram_of(sites, domain);

    HessianBlocks blocks = CapacityHessian<FieldType>(field).assemble(*diagram, weights);

    for (int trial = 0; trial < 3; ++trial) {
        std::vector<Vector3> direction(count);
        std::vector<Vector2> forward = sites;
        std::vector<Vector2> backward = sites;

        for (size_t k = 0; k < count; ++k) {
            direction[k] = Vector3(
                std::sin(1.3 * static_cast<double>(k) + trial),
                std::cos(2.1 * static_cast<double>(k) - trial),
                0.0
            );
            forward[k] += DISPLACEMENT * Vector2(direction[k].x(), direction[k].y());
            backward[k] -= DISPLACEMENT * Vector2(direction[k].x(), direction[k].y());
        }

        std::vector<Vector3> analytic = blocks.multiply(direction);
        std::vector<Vector3> ahead = weighted_capacity_gradient(field, forward, domain, weights);
        std::vector<Vector3> behind = weighted_capacity_gradient(field, backward, domain, weights);

        double scale = 0.0;

        for (size_t k = 0; k < count; ++k) {
            scale = std::max(scale, ((ahead[k] - behind[k]) / (2.0 * DISPLACEMENT)).norm());
        }

        for (size_t k = 0; k < count; ++k) {
            Vector3 difference = (ahead[k] - behind[k]) / (2.0 * DISPLACEMENT);

            EXPECT_LT((analytic[k] - difference).norm(), tolerance * scale) <<
                "site " << k << " trial " << trial;
        }
    }
}

template<fields::flat::Field FieldType>
void expect_symmetric(const FieldType& field, const Domain& domain) {
    std::vector<Vector2> sites = flat_scatter(13, domain.width, domain.height);
    auto diagram = diagram_of(sites, domain);
    std::vector<double> weights = arbitrary_weights(13);

    HessianBlocks blocks = CapacityHessian<FieldType>(field).assemble(*diagram, weights);

    std::vector<Vector3> first(sites.size());
    std::vector<Vector3> second(sites.size());

    for (size_t k = 0; k < sites.size(); ++k) {
        first[k] = Vector3(std::sin(1.1 * static_cast<double>(k)), std::cos(0.6 * static_cast<double>(k)), 0.0);
        second[k] = Vector3(0.2, std::sin(2.3 * static_cast<double>(k)), 0.0);
    }

    std::vector<Vector3> applied_first = blocks.multiply(first);
    std::vector<Vector3> applied_second = blocks.multiply(second);

    double forward = 0.0;
    double backward = 0.0;
    double scale = 0.0;

    for (size_t k = 0; k < sites.size(); ++k) {
        forward += second[k].dot(applied_first[k]);
        backward += first[k].dot(applied_second[k]);
        scale = std::max(scale, applied_first[k].norm());
    }

    EXPECT_NEAR(forward, backward, 1e-10 * std::max(scale, 1.0));
}

} // namespace

// The reference every curvature model is held to: the derivative of the
// exact weighted-capacity gradient. On a constant density the line sweep
// vanishes and only the endpoint, separation and direct terms remain.
TEST(FlatCapacityHessianTest, MatchesFiniteDifferencesOnAConstantField) {
    Domain domain = Domain::torus(1.0, 1.0);
    expect_matches_finite_differences(PolynomialField::constant(1.0, domain), 12, domain, 1e-4);
}

TEST(FlatCapacityHessianTest, MatchesFiniteDifferencesOnAnElongatedTorus) {
    Domain domain = Domain::torus(3.0, 0.5);
    expect_matches_finite_differences(PolynomialField::constant(1.0, domain), 8, domain, 1e-4);
}

TEST(FlatCapacityHessianTest, EXPENSIVE_MatchesFiniteDifferencesOnNoise) {
    REQUIRE_EXPENSIVE();

    Domain domain = Domain::torus(1.0, 1.0);
    NoiseField noise(1.0, 1.0, Interval(0.2, 1.0));
    auto field = PiecewisePolynomialField::sample(domain, 8, 8, noise);

    expect_matches_finite_differences(field, 14, domain, 1e-3);
}

TEST(FlatCapacityHessianTest, IsSymmetric) {
    Domain domain = Domain::torus(1.0, 1.0);
    expect_symmetric(PolynomialField::constant(1.0, domain), domain);
}

// On the plane the outer cells end at the walls, so bisector endpoints on
// a wall slide along it as the sites move instead of following a third
// site: the wall-pinned endpoint velocity is what these exercise.
TEST(PlaneCapacityHessianTest, MatchesFiniteDifferencesOnAConstantField) {
    Domain domain = Domain::plane(1.0, 1.0);
    expect_matches_finite_differences(PolynomialField::constant(1.0, domain), 12, domain, 1e-4);
}

// A gradient adds the line sweep's density terms, felt at the walls too.
TEST(PlaneCapacityHessianTest, MatchesFiniteDifferencesOnALinearField) {
    Domain domain = Domain::plane(2.0, 1.0);
    expect_matches_finite_differences(PolynomialField::linear(1.0, Vector2(1.5, 0.0), domain), 12, domain, 1e-4);
}

TEST(PlaneCapacityHessianTest, EXPENSIVE_MatchesFiniteDifferencesOnNoise) {
    REQUIRE_EXPENSIVE();

    Domain domain = Domain::plane(1.0, 1.0);
    NoiseField noise(1.0, 1.0, Interval(0.2, 1.0));
    auto field = PiecewisePolynomialField::sample(domain, 8, 8, noise);

    expect_matches_finite_differences(field, 14, domain, 1e-3);
}

// The cylinder mixes the two closures: seam neighbours a period apart on
// one axis, wall-pinned endpoints on the other, and a gradient up the
// walled height.
TEST(CylinderCapacityHessianTest, MatchesFiniteDifferencesOnAGradient) {
    Domain domain = Domain::cylinder(2.0, 1.0);
    expect_matches_finite_differences(PolynomialField::linear(1.0, Vector2(0.0, 2.0), domain), 12, domain, 1e-4);
}

TEST(PlaneCapacityHessianTest, IsSymmetric) {
    Domain domain = Domain::plane(1.0, 1.0);
    expect_symmetric(PolynomialField::linear(1.0, Vector2(2.0, 0.0), domain), domain);
}
