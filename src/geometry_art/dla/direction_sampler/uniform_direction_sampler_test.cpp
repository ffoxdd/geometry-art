#include "gtest/gtest.h"
#include "direction_sampler.hpp"
#include "uniform_direction_sampler.hpp"
#include "../../math/interval_sampler/uniform_interval_sampler.hpp"
#include "../../testing/assertions/statistical.hpp"
#include "../../testing/macros.hpp"
#include "../../testing/mocks/interval_sampler.hpp"
#include "../../types.hpp"
#include <vector>

using namespace geometry_art::dla;
using geometry_art::UniformIntervalSampler;
using geometry_art::Vector3;
using geometry_art::VectorS2;
using geometry_art::testing::compute_coordinate_statistics;
using geometry_art::testing::mocks::MockIntervalSampler;

TEST(UniformDirectionSamplerTest, HeightThenAzimuthDrawsProduceTheDirection) {
    UniformDirectionSampler<MockIntervalSampler> sampler(MockIntervalSampler({0.5, 0.0}));

    Vector3 direction = sampler.sample();

    EXPECT_NEAR(direction.x(), 1.0, 1e-12);
    EXPECT_NEAR(direction.y(), 0.0, 1e-12);
    EXPECT_NEAR(direction.z(), 0.0, 1e-12);
}

TEST(UniformDirectionSamplerTest, FullHeightDrawProducesThePole) {
    UniformDirectionSampler<MockIntervalSampler> sampler(MockIntervalSampler({1.0, 0.0}));

    Vector3 direction = sampler.sample();

    EXPECT_NEAR((direction - Vector3(0.0, 0.0, 1.0)).norm(), 0.0, 1e-12);
}

TEST(UniformDirectionSamplerTest, SampleIsUnitLength) {
    UniformDirectionSampler<> sampler{UniformIntervalSampler(7)};

    EXPECT_NEAR(sampler.sample().norm(), 1.0, 1e-12);
}

TEST(UniformDirectionSamplerTest, EXPENSIVE_MeanIsTheCenter) {
    REQUIRE_EXPENSIVE();

    UniformDirectionSampler<> sampler;
    constexpr size_t sample_count = 50000;

    std::vector<VectorS2> directions;
    directions.reserve(sample_count);

    for (size_t i = 0; i < sample_count; ++i) {
        directions.push_back(sampler.sample());
    }

    auto metrics = compute_coordinate_statistics(directions);

    EXPECT_NEAR(metrics.x.mean, 0.0, 0.02);
    EXPECT_NEAR(metrics.y.mean, 0.0, 0.02);
    EXPECT_NEAR(metrics.z.mean, 0.0, 0.02);
}
