#include "gtest/gtest.h"
#include "particle_index.hpp"
#include "rtree_particle_index.hpp"
#include "../../types.hpp"

using namespace geometry_art::dla;
using geometry_art::Vector3;

TEST(RTreeParticleIndexTest, NearestOfASinglePoint) {
    RTreeParticleIndex index;
    index.insert(Vector3(0.0, 0.0, 0.0), 0);

    NearestParticle nearest = index.nearest(Vector3(0.0, 0.0, 3.0));

    EXPECT_EQ(nearest.index, 0u);
    EXPECT_DOUBLE_EQ(nearest.distance, 3.0);
}

TEST(RTreeParticleIndexTest, NearestPicksTheCloserOfTwo) {
    RTreeParticleIndex index;
    index.insert(Vector3(0.0, 0.0, 0.0), 0);
    index.insert(Vector3(0.0, 0.0, 10.0), 1);

    NearestParticle nearest = index.nearest(Vector3(0.0, 0.0, 7.0));

    EXPECT_EQ(nearest.index, 1u);
    EXPECT_DOUBLE_EQ(nearest.distance, 3.0);
}

TEST(RTreeParticleIndexTest, DistanceIsEuclidean) {
    RTreeParticleIndex index;
    index.insert(Vector3(1.0, 2.0, 2.0), 5);

    NearestParticle nearest = index.nearest(Vector3(0.0, 0.0, 0.0));

    EXPECT_EQ(nearest.index, 5u);
    EXPECT_DOUBLE_EQ(nearest.distance, 3.0);
}
