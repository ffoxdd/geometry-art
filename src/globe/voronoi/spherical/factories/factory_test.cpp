#include "factory.hpp"
#include "../../../fields/spherical/polynomial_field.hpp"
#include "../../../testing/macros.hpp"
#include <gtest/gtest.h>

using namespace globe;
using namespace globe::voronoi::spherical;

namespace {

CapacityConstrainedParameters quick_parameters() {
    CapacityConstrainedParameters parameters;
    parameters.max_outer_iterations = 2;
    parameters.max_inner_iterations = 10;
    return parameters;
}

}

constexpr unsigned int TEST_SEED = 20250822;

TEST(FactoryTest, BuildWithConstantDensityPreservesPointCount) {
    Factory factory(6, "constant", 1, "lloyd", 0, quick_parameters(), TEST_SEED, noop_callback());

    auto sphere = factory.build();

    EXPECT_EQ(sphere->size(), 6u);
}

TEST(FactoryTest, BuildWithNewtonWarmStartPreservesPointCount) {
    Factory factory(6, "constant", 0, "newton", 3, quick_parameters(), TEST_SEED, noop_callback());

    auto sphere = factory.build();

    EXPECT_EQ(sphere->size(), 6u);
}

TEST(FactoryTest, EXPENSIVE_BuildWithQuadraticDensityEqualizesCapacities) {
    REQUIRE_EXPENSIVE();

    CapacityConstrainedParameters parameters;
    parameters.relative_capacity_tolerance = 1e-6;
    Factory factory(40, "quadratic", 5, "lloyd", 0, parameters, TEST_SEED, noop_callback());

    auto sphere = factory.build();

    Eigen::Matrix3d quadratic = Eigen::Matrix3d::Zero();
    quadratic(2, 2) = -0.9;
    auto field = fields::spherical::PolynomialField::quadratic(1.0, Vector3::Zero(), quadratic);
    double target = field.total_mass() / sphere->size();
    double sum = 0.0;
    for (const auto& cell : sphere->cells()) {
        double error = field.integrals(cell).mass - target;
        sum += error * error;
    }

    EXPECT_LT(std::sqrt(sum / sphere->size()) / target, 1e-6);
}

TEST(FactoryTest, EXPENSIVE_BuildWithNoiseDensity) {
    REQUIRE_EXPENSIVE();

    Factory factory(20, "noise", 2, "lloyd", 0, quick_parameters(), TEST_SEED, noop_callback());

    auto sphere = factory.build();

    EXPECT_EQ(sphere->size(), 20u);
}

TEST(FactoryTest, EXPENSIVE_BuildWithFittedNoiseDensity) {
    REQUIRE_EXPENSIVE();

    Factory factory(20, "noise-fit", 2, "lloyd", 0, quick_parameters(), TEST_SEED, noop_callback());

    auto sphere = factory.build();

    EXPECT_EQ(sphere->size(), 20u);
}
