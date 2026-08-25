#include "gtest/gtest.h"
#include "uniform_circular_interval_sampler.hpp"
#include "../circular_interval.hpp"
#include "../../geometry/spherical/helpers.hpp"
#include "../../testing/assertions/geometric.hpp"
#include "../../testing/assertions/statistical.hpp"
#include "../../testing/macros.hpp"

using namespace geometry_art;
using geometry_art::testing::compute_statistics;
using geometry_art::testing::expect_mean;
using geometry_art::testing::expect_variance;
using geometry_art::testing::expect_range_coverage;
using geometry_art::testing::uniform_distribution_mean;
using geometry_art::testing::uniform_distribution_variance;

using ThetaInterval = CircularInterval<TWO_PI>;

TEST(UniformCircularIntervalSamplerTest, SampleFallsWithinSimpleInterval) {
    UniformCircularIntervalSampler sampler(42);
    ThetaInterval interval(M_PI / 4, M_PI / 2);

    double value = sampler.sample(interval);

    EXPECT_GE(value, interval.start());
    EXPECT_LT(value, interval.end());
}

TEST(UniformCircularIntervalSamplerTest, SampleFromWrappedIntervalIsContained) {
    UniformCircularIntervalSampler sampler(42);
    ThetaInterval interval(5.5, 1.5);

    double value = sampler.sample(interval);

    EXPECT_TRUE(interval.contains(value));
}

TEST(UniformCircularIntervalSamplerTest, SampleFromFullIntervalFallsWithinPeriod) {
    UniformCircularIntervalSampler sampler(42);
    ThetaInterval interval = ThetaInterval::full();

    double value = sampler.sample(interval);

    EXPECT_GE(value, 0.0);
    EXPECT_LT(value, TWO_PI);
}

TEST(UniformCircularIntervalSamplerTest, SeededSamplerIsDeterministic) {
    UniformCircularIntervalSampler sampler1(123);
    UniformCircularIntervalSampler sampler2(123);
    ThetaInterval interval(0.0, M_PI);

    double value1 = sampler1.sample(interval);
    double value2 = sampler2.sample(interval);

    EXPECT_DOUBLE_EQ(value1, value2);
}

TEST(UniformCircularIntervalSamplerTest, EXPENSIVE_AllSamplesWithinBounds) {
    REQUIRE_EXPENSIVE();

    UniformCircularIntervalSampler sampler;
    ThetaInterval interval(M_PI / 4, M_PI / 2);
    constexpr size_t sample_count = 10000;

    for (size_t i = 0; i < sample_count; ++i) {
        double value = sampler.sample(interval);
        EXPECT_TRUE(interval.contains(value));
    }
}

TEST(UniformCircularIntervalSamplerTest, EXPENSIVE_UniformDistribution) {
    REQUIRE_EXPENSIVE();

    UniformCircularIntervalSampler sampler;
    ThetaInterval interval(1.0, 2.0);
    constexpr size_t sample_count = 50000;

    auto stats = compute_statistics(
        [&]() { return sampler.sample(interval); },
        sample_count
    );

    double expected_mean = interval.start() + interval.measure() / 2.0;
    double expected_variance = (interval.measure() * interval.measure()) / 12.0;

    expect_mean(stats, expected_mean, 0.05);
    expect_variance(stats, expected_variance, 0.02);
}

TEST(UniformCircularIntervalSamplerTest, EXPENSIVE_WrappedIntervalContainment) {
    REQUIRE_EXPENSIVE();

    UniformCircularIntervalSampler sampler;
    ThetaInterval interval(5.5, 1.5);
    constexpr size_t sample_count = 10000;

    for (size_t i = 0; i < sample_count; ++i) {
        double value = sampler.sample(interval);
        EXPECT_TRUE(interval.contains(value));
    }
}
