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
    Factory factory(6, "constant", 1, "lloyd", 0, quick_parameters(), TEST_SEED, noop_callback(), {}, std::chrono::milliseconds(0));

    auto sphere = factory.build();

    EXPECT_EQ(sphere->size(), 6u);
}

TEST(FactoryTest, BuildWithNewtonWarmStartPreservesPointCount) {
    Factory factory(6, "constant", 0, "newton", 3, quick_parameters(), TEST_SEED, noop_callback(), {}, std::chrono::milliseconds(0));

    auto sphere = factory.build();

    EXPECT_EQ(sphere->size(), 6u);
}

TEST(FactoryTest, EXPENSIVE_BuildWithQuadraticDensityEqualizesCapacities) {
    REQUIRE_EXPENSIVE();

    CapacityConstrainedParameters parameters;
    parameters.relative_capacity_tolerance = 1e-6;
    Factory factory(40, "quadratic", 5, "lloyd", 0, parameters, TEST_SEED, noop_callback(), {}, std::chrono::milliseconds(0));

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

    Factory factory(20, "noise", 2, "lloyd", 0, quick_parameters(), TEST_SEED, noop_callback(), {}, std::chrono::milliseconds(0));

    auto sphere = factory.build();

    EXPECT_EQ(sphere->size(), 20u);
}

TEST(FactoryTest, EXPENSIVE_BuildWithFittedNoiseDensity) {
    REQUIRE_EXPENSIVE();

    Factory factory(20, "noise-fit", 2, "lloyd", 0, quick_parameters(), TEST_SEED, noop_callback(), {}, std::chrono::milliseconds(0));

    auto sphere = factory.build();

    EXPECT_EQ(sphere->size(), 20u);
}

// The mesh follows the site count by the bandwidth rule: the coarsest
// icosphere whose edges fit inside a cell of angular radius 2 / sqrt(N).
TEST(FactoryTest, SmoothNoiseSubdivisionsFollowTheCellScale) {
    EXPECT_EQ(Factory::smooth_noise_subdivisions(10), 1);
    EXPECT_EQ(Factory::smooth_noise_subdivisions(50), 2);
    EXPECT_EQ(Factory::smooth_noise_subdivisions(200), 3);
    EXPECT_EQ(Factory::smooth_noise_subdivisions(500), 4);
    EXPECT_EQ(Factory::smooth_noise_subdivisions(1000), 5);
}

TEST(FactoryTest, SmoothNoiseSubdivisionsAreMonotoneAndClamped) {
    int previous = 0;

    for (int count : {2, 10, 100, 1000, 10000, 100000, 1000000}) {
        int level = Factory::smooth_noise_subdivisions(count);
        EXPECT_GE(level, previous);
        EXPECT_GE(level, 1);
        EXPECT_LE(level, 6);
        previous = level;
    }
}
