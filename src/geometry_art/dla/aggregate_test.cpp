#include "gtest/gtest.h"
#include "aggregate.hpp"
#include "particle.hpp"
#include "../types.hpp"
#include <optional>

using namespace geometry_art::dla;
using geometry_art::Vector3;

namespace {

Aggregate<> seeded_aggregate() {
    return Aggregate<>(1.0, Particle{Vector3::Zero(), std::nullopt});
}

} // namespace

TEST(AggregateTest, StartsWithTheSeedParticle) {
    Aggregate<> aggregate = seeded_aggregate();

    EXPECT_EQ(aggregate.size(), 1u);
    EXPECT_FALSE(aggregate.particles()[0].parent.has_value());
    EXPECT_DOUBLE_EQ(aggregate.particle_radius(), 1.0);
    EXPECT_DOUBLE_EQ(aggregate.contact_distance(), 2.0);
    EXPECT_DOUBLE_EQ(aggregate.reach(), 2.0);
}

TEST(AggregateTest, ClearanceIsTheDistanceToTheContactSphere) {
    Aggregate<> aggregate = seeded_aggregate();

    EXPECT_DOUBLE_EQ(aggregate.clearance(Vector3(0.0, 0.0, 5.0)), 3.0);
}

TEST(AggregateTest, ClearanceUsesTheNearestParticle) {
    Aggregate<> aggregate = seeded_aggregate();
    aggregate.freeze(Particle{Vector3(0.0, 0.0, 4.0), 0});

    EXPECT_DOUBLE_EQ(aggregate.clearance(Vector3(0.0, 0.0, 5.5)), -0.5);
}

TEST(AggregateTest, NearestParticleReportsIndexAndDistance) {
    Aggregate<> aggregate = seeded_aggregate();
    aggregate.freeze(Particle{Vector3(0.0, 0.0, 4.0), 0});

    NearestParticle nearest = aggregate.nearest_particle(Vector3(0.0, 0.0, 5.5));

    EXPECT_EQ(nearest.index, 1u);
    EXPECT_DOUBLE_EQ(nearest.distance, 1.5);
}

TEST(AggregateTest, FreezeExtendsTheReach) {
    Aggregate<> aggregate = seeded_aggregate();
    aggregate.freeze(Particle{Vector3(0.0, 0.0, 4.0), 0});

    EXPECT_DOUBLE_EQ(aggregate.reach(), 6.0);
}
