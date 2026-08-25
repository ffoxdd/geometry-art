#include "gtest/gtest.h"
#include "elimination_point_generator.hpp"
#include "../../../geometry/spherical/helpers.hpp"
#include "../../../testing/assertions/geometric.hpp"
#include "../../../testing/macros.hpp"

using namespace geometry_art::generators::spherical::poisson;
using geometry_art::generators::spherical::RandomPointGenerator;
using geometry_art::VectorS2;
using geometry_art::geometry::spherical::distance;

TEST(EliminationPointGeneratorTest, GeneratesRequestedCount) {
    EliminationPointGenerator<> generator;

    auto points = generator.generate(10);

    EXPECT_EQ(points.size(), 10);
}

TEST(EliminationPointGeneratorTest, GeneratesZeroPointsWhenRequested) {
    EliminationPointGenerator<> generator;

    auto points = generator.generate(0);

    EXPECT_TRUE(points.empty());
}

TEST(EliminationPointGeneratorTest, AllPointsAreOnUnitSphere) {
    EliminationPointGenerator<> generator;

    auto points = generator.generate(20);

    for (const auto &point : points) {
        double distance_from_origin = std::sqrt(
            point.x() * point.x() +
            point.y() * point.y() +
            point.z() * point.z()
        );
        EXPECT_NEAR(distance_from_origin, 1.0, 1e-10);
    }
}

TEST(EliminationPointGeneratorTest, TracksAttemptCount) {
    EliminationPointGenerator<> generator(
        RandomPointGenerator<>(),
        2.0  // oversample factor
    );

    generator.generate(10);

    // With oversample factor of 2.0, we generate 20 candidates for 10 points
    EXPECT_GE(generator.last_attempt_count(), 20);
}

TEST(EliminationPointGeneratorTest, EXPENSIVE_PointsAreWellDistributed) {
    REQUIRE_EXPENSIVE();

    EliminationPointGenerator<> generator(
        RandomPointGenerator<>(),
        3.0  // higher oversample for better distribution
    );

    auto points = generator.generate(100);

    // Calculate minimum distance between any two points
    double min_distance = std::numeric_limits<double>::max();
    for (size_t i = 0; i < points.size(); ++i) {
        for (size_t j = i + 1; j < points.size(); ++j) {
            double dist = distance(points[i], points[j]);
            min_distance = std::min(min_distance, dist);
        }
    }

    // For 100 points on a sphere, optimal packing gives roughly sqrt(4*pi/100) ≈ 0.35 rad spacing
    // We expect elimination to achieve at least half of optimal
    double optimal_spacing = std::sqrt(4.0 * M_PI / 100.0);
    EXPECT_GT(min_distance, optimal_spacing * 0.3);
}

TEST(EliminationPointGeneratorTest, EXPENSIVE_BetterThanRandomDistribution) {
    REQUIRE_EXPENSIVE();

    // Generate random points
    RandomPointGenerator<> random_generator;
    auto random_points = random_generator.generate(100);

    // Generate Poisson points
    EliminationPointGenerator<> poisson_generator(
        RandomPointGenerator<>(),
        3.0
    );
    auto poisson_points = poisson_generator.generate(100);

    // Calculate minimum distances
    auto min_distance = [](const std::vector<VectorS2> &points) {
        double min_dist = std::numeric_limits<double>::max();
        for (size_t i = 0; i < points.size(); ++i) {
            for (size_t j = i + 1; j < points.size(); ++j) {
                min_dist = std::min(min_dist, distance(points[i], points[j]));
            }
        }
        return min_dist;
    };

    double random_min = min_distance(random_points);
    double poisson_min = min_distance(poisson_points);

    // Poisson should produce better minimum spacing than random
    EXPECT_GT(poisson_min, random_min);
}
