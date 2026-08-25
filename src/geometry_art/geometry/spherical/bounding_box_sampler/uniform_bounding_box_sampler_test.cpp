#include "gtest/gtest.h"
#include "bounding_box_sampler.hpp"
#include "uniform_bounding_box_sampler.hpp"
#include "../bounding_box.hpp"
#include "../../../math/interval.hpp"
#include "../../../testing/assertions/geometric.hpp"
#include "../../../testing/assertions/statistical.hpp"
#include "../../../testing/macros.hpp"
#include "../../../testing/mocks/interval_sampler.hpp"
#include "../../../testing/mocks/circular_interval_sampler.hpp"
#include <cmath>

using namespace geometry_art;
using geometry_art::testing::compute_statistics;
using geometry_art::testing::expect_mean;
using geometry_art::testing::expect_variance;
using geometry_art::testing::expect_range_coverage;
using geometry_art::testing::uniform_distribution_mean;
using geometry_art::testing::uniform_distribution_variance;
using geometry_art::testing::mocks::MockIntervalSampler;
using geometry_art::testing::mocks::MockCircularIntervalSampler;

TEST(UniformSphericalBoundingBoxSamplerTest, ProducesPointOnUnitSphere) {
    MockIntervalSampler z_sampler({0.5});
    MockCircularIntervalSampler theta_sampler({0.25});
    UniformSphericalBoundingBoxSampler sampler(z_sampler, theta_sampler);
    SphericalBoundingBox box{ThetaInterval::full(), Interval(-1.0, 1.0)};

    VectorS2 sample = sampler.sample(box);

    EXPECT_NEAR(sample.norm(), 1.0, 1e-9);
}

TEST(UniformSphericalBoundingBoxSamplerTest, UsesZFromIntervalSampler) {
    MockIntervalSampler z_sampler({0.75});
    MockCircularIntervalSampler theta_sampler({0.0});
    UniformSphericalBoundingBoxSampler sampler(z_sampler, theta_sampler);
    SphericalBoundingBox box{ThetaInterval::full(), Interval(0.0, 1.0)};

    VectorS2 sample = sampler.sample(box);

    EXPECT_NEAR(sample.z(), 0.75, 1e-9);
}

TEST(UniformSphericalBoundingBoxSamplerTest, UsesThetaFromCircularIntervalSampler) {
    MockIntervalSampler z_sampler({0.5});
    MockCircularIntervalSampler theta_sampler({0.5});
    UniformSphericalBoundingBoxSampler sampler(z_sampler, theta_sampler);
    SphericalBoundingBox box{ThetaInterval(0.0, M_PI), Interval(-1.0, 1.0)};

    VectorS2 sample = sampler.sample(box);
    double theta = std::atan2(sample.y(), sample.x());

    EXPECT_NEAR(theta, M_PI / 2.0, 1e-9);
}

TEST(UniformSphericalBoundingBoxSamplerTest, HandlesWrappedThetaInterval) {
    MockIntervalSampler z_sampler({0.5});
    MockCircularIntervalSampler theta_sampler({0.0});
    UniformSphericalBoundingBoxSampler sampler(z_sampler, theta_sampler);
    SphericalBoundingBox box{ThetaInterval(5.5, 1.5), Interval(-1.0, 1.0)};

    VectorS2 sample = sampler.sample(box);

    EXPECT_TRUE(box.contains(sample));
}

TEST(UniformSphericalBoundingBoxSamplerTest, EXPENSIVE_AllPointsOnUnitSphere) {
    REQUIRE_EXPENSIVE();

    UniformSphericalBoundingBoxSampler<> sampler;
    SphericalBoundingBox box{ThetaInterval::full(), Interval(0.0, 0.5)};
    constexpr size_t sample_count = 10000;

    for (size_t i = 0; i < sample_count; ++i) {
        VectorS2 sample = sampler.sample(box);

        EXPECT_NEAR(sample.norm(), 1.0, 1e-9);
        EXPECT_TRUE(box.contains(sample));
    }
}

TEST(UniformSphericalBoundingBoxSamplerTest, EXPENSIVE_UniformZDistribution) {
    REQUIRE_EXPENSIVE();

    UniformSphericalBoundingBoxSampler<> sampler;
    Interval z_interval(-0.3, 0.7);
    SphericalBoundingBox box{ThetaInterval::full(), z_interval};
    constexpr size_t sample_count = 50000;

    auto stats = compute_statistics(
        [&]() { return sampler.sample(box).z(); },
        sample_count
    );

    expect_mean(stats, uniform_distribution_mean(z_interval.low(), z_interval.high()), 0.02);
    expect_variance(stats, uniform_distribution_variance(z_interval.low(), z_interval.high()), 0.01);
    expect_range_coverage(stats, z_interval, 0.95);
}

TEST(UniformSphericalBoundingBoxSamplerTest, EXPENSIVE_WrappedIntervalContainment) {
    REQUIRE_EXPENSIVE();

    UniformSphericalBoundingBoxSampler<> sampler;
    SphericalBoundingBox box{ThetaInterval(5.5, 1.58), Interval(-0.2, 0.2)};
    constexpr size_t sample_count = 10000;

    for (size_t i = 0; i < sample_count; ++i) {
        VectorS2 sample = sampler.sample(box);

        EXPECT_NEAR(sample.norm(), 1.0, 1e-9);
        EXPECT_TRUE(box.contains(sample));
    }
}

TEST(UniformSphericalBoundingBoxSamplerTest, EXPENSIVE_CrossingZeroContainment) {
    REQUIRE_EXPENSIVE();

    UniformSphericalBoundingBoxSampler<> sampler;
    SphericalBoundingBox box{ThetaInterval(2 * M_PI - 0.5, 1.0), Interval(0.0, 0.5)};
    constexpr size_t sample_count = 10000;

    for (size_t i = 0; i < sample_count; ++i) {
        VectorS2 sample = sampler.sample(box);

        EXPECT_TRUE(box.contains(sample));
    }
}
