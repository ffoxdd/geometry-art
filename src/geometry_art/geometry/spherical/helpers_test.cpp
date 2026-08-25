#include <gtest/gtest.h>
#include "helpers.hpp"
#include <cmath>

using namespace geometry_art;
using geometry::spherical::distance;

TEST(DistanceTest, ComputesRightAngle) {
    VectorS2 a(1, 0, 0);
    VectorS2 b(0, 1, 0);

    double angle = distance(a, b);

    EXPECT_DOUBLE_EQ(angle, M_PI_2);
}

TEST(DistanceTest, ComputesOppositeVectors) {
    VectorS2 a(1, 0, 0);
    VectorS2 b(-1, 0, 0);

    double angle = distance(a, b);

    EXPECT_DOUBLE_EQ(angle, M_PI);
}

TEST(DistanceTest, ComputesSameDirection) {
    VectorS2 a(1, 0, 0);
    VectorS2 b(1, 0, 0);

    double angle = distance(a, b);

    EXPECT_NEAR(angle, 0.0, 1e-10);
}
