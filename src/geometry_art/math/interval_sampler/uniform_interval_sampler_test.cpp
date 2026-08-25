#include "gtest/gtest.h"
#include "interval_sampler.hpp"
#include "uniform_interval_sampler.hpp"
#include "../interval.hpp"
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

TEST(UniformIntervalSamplerTest, SampleFallsWithinInterval) {
    UniformIntervalSampler sampler(42);
    Interval interval(2.5, 5.0);

    double value = sampler.sample(interval);

    EXPECT_GE(value, interval.low());
    EXPECT_LT(value, interval.high());
}

TEST(UniformIntervalSamplerTest, SampleFromDegenerateIntervalReturnsLow) {
    UniformIntervalSampler sampler(42);
    Interval interval(3.0, 3.0);

    double value = sampler.sample(interval);

    EXPECT_DOUBLE_EQ(value, 3.0);
}

TEST(UniformIntervalSamplerTest, SeededSamplerIsDeterministic) {
    UniformIntervalSampler sampler1(123);
    UniformIntervalSampler sampler2(123);
    Interval interval(0.0, 1.0);

    double value1 = sampler1.sample(interval);
    double value2 = sampler2.sample(interval);

    EXPECT_DOUBLE_EQ(value1, value2);
}

TEST(UniformIntervalSamplerTest, EXPENSIVE_AllSamplesWithinBounds) {
    REQUIRE_EXPENSIVE();

    UniformIntervalSampler sampler;
    Interval interval(2.5, 5.0);
    constexpr size_t sample_count = 10000;

    for (size_t i = 0; i < sample_count; ++i) {
        double value = sampler.sample(interval);
        EXPECT_GE(value, interval.low());
        EXPECT_LT(value, interval.high());
    }
}

TEST(UniformIntervalSamplerTest, EXPENSIVE_UniformDistribution) {
    REQUIRE_EXPENSIVE();

    UniformIntervalSampler sampler;
    Interval interval(-3.0, 7.0);
    constexpr size_t sample_count = 50000;

    auto stats = compute_statistics(
        [&]() { return sampler.sample(interval); },
        sample_count
    );

    expect_mean(stats, uniform_distribution_mean(interval.low(), interval.high()), 0.1);
    expect_variance(stats, uniform_distribution_variance(interval.low(), interval.high()), 0.5);
    expect_range_coverage(stats, interval, 0.95);
}
