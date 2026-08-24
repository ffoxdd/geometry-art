#include "capacity_hessian.hpp"
#include "capacity_constrained_lagrangian.hpp"
#include "../core/torus.hpp"
#include "../../capacity_jacobian.hpp"
#include "../../hessian_blocks.hpp"
#include "../../../fields/flat/constant_field.hpp"
#include "../../../fields/flat/noise_field.hpp"
#include "../../../fields/flat/piecewise_polynomial_field.hpp"
#include "../../../math/interval.hpp"
#include "../../../testing/macros.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

using namespace globe;
using namespace globe::voronoi;
using namespace globe::voronoi::flat;
using fields::flat::ConstantField;
using fields::flat::NoiseField;
using fields::flat::PiecewisePolynomialField;
using globe::math::Interval;

namespace {

constexpr double DISPLACEMENT = 1e-6;

// A two-dimensional low-discrepancy scatter. Rows of sites make nearly
// striped cells whose capacities are close to linear in the sites, which
// starves a curvature test of its subject.
std::vector<Vector2> scattered_sites(size_t count, double width, double height) {
    std::vector<Vector2> sites;

    for (size_t k = 0; k < count; ++k) {
        double x = std::fmod(0.13 + 0.7548776662466927 * static_cast<double>(k), 1.0) * width;
        double y = std::fmod(0.41 + 0.5698402909980532 * static_cast<double>(k), 1.0) * height;
        sites.emplace_back(x, y);
    }

    return sites;
}

std::unique_ptr<Torus> torus_of(const std::vector<Vector2>& sites, double width, double height) {
    auto torus = std::make_unique<Torus>(width, height);

    for (const Vector2& site : sites) {
        torus->insert(site);
    }

    return torus;
}

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
    double width,
    double height,
    const std::vector<double>& weights
) {
    CapacityConstrainedLagrangian<FieldType> lagrangian(field, 0.0);
    auto torus = torus_of(sites, width, height);

    return CapacityJacobian(lagrangian.diagram_state(*torus)).transpose_apply(weights);
}

template<fields::flat::Field FieldType>
void expect_matches_finite_differences(
    const FieldType& field,
    size_t count,
    double width,
    double height,
    double tolerance
) {
    std::vector<Vector2> sites = scattered_sites(count, width, height);
    std::vector<double> weights = arbitrary_weights(count);
    auto torus = torus_of(sites, width, height);

    HessianBlocks blocks = CapacityHessian<FieldType>(field).assemble(*torus, weights);

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
        std::vector<Vector3> ahead = weighted_capacity_gradient(field, forward, width, height, weights);
        std::vector<Vector3> behind = weighted_capacity_gradient(field, backward, width, height, weights);

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

} // namespace

// The reference every curvature model is held to: the derivative of the
// exact weighted-capacity gradient. On a constant density the line sweep
// vanishes and only the endpoint, separation and direct terms remain.
TEST(FlatCapacityHessianTest, MatchesFiniteDifferencesOnAConstantField) {
    expect_matches_finite_differences(ConstantField(1.0, 1.0, 1.0), 12, 1.0, 1.0, 1e-4);
}

TEST(FlatCapacityHessianTest, MatchesFiniteDifferencesOnAnElongatedTorus) {
    expect_matches_finite_differences(ConstantField(1.0, 3.0, 0.5), 8, 3.0, 0.5, 1e-4);
}

TEST(FlatCapacityHessianTest, EXPENSIVE_MatchesFiniteDifferencesOnNoise) {
    REQUIRE_EXPENSIVE();

    NoiseField noise(1.0, 1.0, Interval(0.2, 1.0));
    auto field = PiecewisePolynomialField::sample(1.0, 1.0, 8, 8, noise);

    expect_matches_finite_differences(field, 14, 1.0, 1.0, 1e-3);
}

TEST(FlatCapacityHessianTest, IsSymmetric) {
    ConstantField field(1.0, 1.0, 1.0);
    std::vector<Vector2> sites = scattered_sites(13, 1.0, 1.0);
    auto torus = torus_of(sites, 1.0, 1.0);
    std::vector<double> weights = arbitrary_weights(13);

    HessianBlocks blocks = CapacityHessian<ConstantField>(field).assemble(*torus, weights);

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
