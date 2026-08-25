#include <gtest/gtest.h>
#include "noise_field.hpp"
#include "../../testing/assertions/geometric.hpp"
#include "../../testing/assertions/statistical.hpp"
#include "../../testing/macros.hpp"
#include "../../math/interval.hpp"
#include "../../generators/spherical/random_point_generator.hpp"

using namespace geometry_art::fields::scalar;
using geometry_art::Interval;
using geometry_art::generators::spherical::RandomPointGenerator;
using geometry_art::VectorS2;
using geometry_art::testing::compute_statistics;
using geometry_art::testing::expect_range_coverage;

TEST(NoiseFieldTest, ValueMethodReturnsConsistentResult) {
    VectorS2 location = {0.1, 0.2, 0.3};
    NoiseField noise_field;

    double value1 = noise_field.value(location);
    double value2 = noise_field.value(location);

    EXPECT_DOUBLE_EQ(value1, value2);
}

TEST(NoiseFieldTest, ValueMethodDifferentLocations) {
    VectorS2 location1 = {0.1, 0.2, 0.3};
    VectorS2 location2 = {0.4, 0.5, 0.6};
    NoiseField noise_field;

    double value1 = noise_field.value(location1);
    double value2 = noise_field.value(location2);

    EXPECT_NE(value1, value2);
}

constexpr int SAMPLE_COUNT = 200;
constexpr int TEST_SEED_COUNT = 5000;

std::vector<int> test_seeds() {
    std::vector<int> seeds;
    for (int i = 0; i < TEST_SEED_COUNT; ++i) {
        seeds.push_back(i * 137 + 42);
    }
    return seeds;
}

TEST(NoiseFieldTest, EXPENSIVE_CanConfigureOutputRange) {
    REQUIRE_EXPENSIVE();

    Interval output_range(-0.01, 0.02);
    RandomPointGenerator point_generator;

    for (int seed : test_seeds()) {
        NoiseField noise_field(output_range, seed);

        auto points = point_generator.generate(SAMPLE_COUNT);
        std::vector<double> values;
        values.reserve(SAMPLE_COUNT);
        for (const auto& point : points) {
            values.push_back(noise_field.value(point));
        }

        auto metrics = compute_statistics(values);

        EXPECT_GE(metrics.min_value, output_range.low())
            << "Failed for seed " << seed << " (min_value was " << metrics.min_value << ")";

        EXPECT_LE(metrics.max_value, output_range.high())
            << "Failed for seed " << seed << " (max_value was " << metrics.max_value << ")";
    }
}

TEST(NoiseFieldTest, EXPENSIVE_OutputDistributionUsesFullRange) {
    REQUIRE_EXPENSIVE();

    Interval expected_range(0, 1);
    RandomPointGenerator point_generator;

    for (int seed : test_seeds()) {
        NoiseField noise_field(expected_range, seed);

        auto points = point_generator.generate(SAMPLE_COUNT);
        std::vector<double> values;
        values.reserve(SAMPLE_COUNT);
        for (const auto& point : points) {
            values.push_back(noise_field.value(point));
        }

        auto metrics = compute_statistics(values);

        expect_range_coverage(metrics, expected_range, 0.8);
        EXPECT_LT(metrics.clipping_ratio, 0.1) << "Failed for seed " << seed;
    }
}
