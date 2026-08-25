#include "gtest/gtest.h"
#include "bounding_box_sampler.hpp"
#include "uniform_bounding_box_sampler.hpp"
#include "../bounding_box.hpp"
#include "../../../math/interval.hpp"
#include "../../../testing/assertions/geometric.hpp"
#include "../../../testing/assertions/statistical.hpp"
#include "../../../testing/macros.hpp"
#include "../../../testing/mocks/interval_sampler.hpp"

using namespace geometry_art;
using geometry_art::testing::compute_statistics;
using geometry_art::testing::expect_mean;
using geometry_art::testing::expect_variance;
using geometry_art::testing::expect_range_coverage;
using geometry_art::testing::uniform_distribution_mean;
using geometry_art::testing::uniform_distribution_variance;
using geometry_art::testing::mocks::MockIntervalSampler;

TEST(UniformBoundingBoxSamplerTest, UsesXFromIntervalSampler) {
    MockIntervalSampler sampler({0.25, 0.0, 0.0});
    UniformBoundingBoxSampler bounding_box_sampler(sampler);
    BoundingBox box(Interval(0.0, 4.0), Interval(0.0, 1.0), Interval(0.0, 1.0));

    Vector3 sample = bounding_box_sampler.sample(box);

    EXPECT_NEAR(sample.x(), 1.0, 1e-9);
}

TEST(UniformBoundingBoxSamplerTest, UsesYFromIntervalSampler) {
    MockIntervalSampler sampler({0.0, 0.5, 0.0});
    UniformBoundingBoxSampler bounding_box_sampler(sampler);
    BoundingBox box(Interval(0.0, 1.0), Interval(2.0, 6.0), Interval(0.0, 1.0));

    Vector3 sample = bounding_box_sampler.sample(box);

    EXPECT_NEAR(sample.y(), 4.0, 1e-9);
}

TEST(UniformBoundingBoxSamplerTest, UsesZFromIntervalSampler) {
    MockIntervalSampler sampler({0.0, 0.0, 0.75});
    UniformBoundingBoxSampler bounding_box_sampler(sampler);
    BoundingBox box(Interval(0.0, 1.0), Interval(0.0, 1.0), Interval(-2.0, 2.0));

    Vector3 sample = bounding_box_sampler.sample(box);

    EXPECT_NEAR(sample.z(), 1.0, 1e-9);
}

TEST(UniformBoundingBoxSamplerTest, SampleFallsWithinBoundingBox) {
    UniformBoundingBoxSampler<> sampler;
    BoundingBox box(Interval(1.0, 3.0), Interval(-2.0, 2.0), Interval(0.0, 5.0));

    Vector3 sample = sampler.sample(box);

    EXPECT_GE(sample.x(), box.x_interval().low());
    EXPECT_LT(sample.x(), box.x_interval().high());
    EXPECT_GE(sample.y(), box.y_interval().low());
    EXPECT_LT(sample.y(), box.y_interval().high());
    EXPECT_GE(sample.z(), box.z_interval().low());
    EXPECT_LT(sample.z(), box.z_interval().high());
}

TEST(UniformBoundingBoxSamplerTest, EXPENSIVE_AllSamplesWithinBounds) {
    REQUIRE_EXPENSIVE();

    UniformBoundingBoxSampler<> sampler;
    BoundingBox box(Interval(1.0, 3.0), Interval(-2.0, 2.0), Interval(0.0, 5.0));
    constexpr size_t sample_count = 10000;

    for (size_t i = 0; i < sample_count; ++i) {
        Vector3 sample = sampler.sample(box);

        EXPECT_GE(sample.x(), box.x_interval().low());
        EXPECT_LT(sample.x(), box.x_interval().high());
        EXPECT_GE(sample.y(), box.y_interval().low());
        EXPECT_LT(sample.y(), box.y_interval().high());
        EXPECT_GE(sample.z(), box.z_interval().low());
        EXPECT_LT(sample.z(), box.z_interval().high());
    }
}

TEST(UniformBoundingBoxSamplerTest, EXPENSIVE_UniformXDistribution) {
    REQUIRE_EXPENSIVE();

    UniformBoundingBoxSampler<> sampler;
    Interval x_interval(1.0, 5.0);
    BoundingBox box(x_interval, Interval(0.0, 1.0), Interval(0.0, 1.0));
    constexpr size_t sample_count = 50000;

    auto stats = compute_statistics(
        [&]() { return sampler.sample(box).x(); },
        sample_count
    );

    expect_mean(stats, uniform_distribution_mean(x_interval.low(), x_interval.high()), 0.05);
    expect_variance(stats, uniform_distribution_variance(x_interval.low(), x_interval.high()), 0.1);
    expect_range_coverage(stats, x_interval, 0.95);
}
