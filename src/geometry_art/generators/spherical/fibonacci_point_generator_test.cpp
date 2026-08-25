#include "fibonacci_point_generator.hpp"
#include "../../geometry/spherical/bounding_box.hpp"
#include "../../math/circular_interval.hpp"
#include "../../math/interval.hpp"
#include <gtest/gtest.h>
#include <cmath>

using namespace geometry_art::generators::spherical;
using geometry_art::SphericalBoundingBox;
using geometry_art::CircularInterval;
using geometry_art::Interval;
using geometry_art::TWO_PI;

TEST(FibonacciPointGeneratorTest, GeneratesRequestedCount) {
    FibonacciPointGenerator generator;

    auto points = generator.generate(100);

    EXPECT_EQ(points.size(), 100);
    EXPECT_EQ(generator.last_attempt_count(), 100);
}

TEST(FibonacciPointGeneratorTest, PointsAreOnUnitSphere) {
    FibonacciPointGenerator generator;
    auto points = generator.generate(500);

    for (const auto& point : points) {
        EXPECT_NEAR(point.norm(), 1.0, 1e-10);
    }
}

TEST(FibonacciPointGeneratorTest, PointsAreUniformlyDistributed) {
    FibonacciPointGenerator generator;
    auto points = generator.generate(1000);

    double z_sum = 0.0;
    for (const auto& point : points) {
        z_sum += point.z();
    }
    double z_mean = z_sum / points.size();

    EXPECT_NEAR(z_mean, 0.0, 0.05);
}

TEST(FibonacciPointGeneratorTest, PointsAreDeterministic) {
    FibonacciPointGenerator generator1;
    FibonacciPointGenerator generator2;

    auto points1 = generator1.generate(100);
    auto points2 = generator2.generate(100);

    ASSERT_EQ(points1.size(), points2.size());
    for (size_t i = 0; i < points1.size(); ++i) {
        EXPECT_DOUBLE_EQ(points1[i].x(), points2[i].x());
        EXPECT_DOUBLE_EQ(points1[i].y(), points2[i].y());
        EXPECT_DOUBLE_EQ(points1[i].z(), points2[i].z());
    }
}

TEST(FibonacciPointGeneratorTest, GenerateWithBoundingBoxFiltersPoints) {
    FibonacciPointGenerator generator;

    SphericalBoundingBox northern_hemisphere(
        CircularInterval<TWO_PI>::full(),
        Interval(0.0, 1.0)
    );

    auto points = generator.generate(1000, northern_hemisphere);

    EXPECT_LT(points.size(), 1000);
    EXPECT_GT(points.size(), 400);
    EXPECT_EQ(generator.last_attempt_count(), 1000);

    for (const auto& point : points) {
        EXPECT_GE(point.z(), -1e-10);
    }
}

TEST(FibonacciPointGeneratorTest, EmptyGenerationReturnsEmpty) {
    FibonacciPointGenerator generator;

    auto points = generator.generate(0);

    EXPECT_TRUE(points.empty());
    EXPECT_EQ(generator.last_attempt_count(), 0);
}
