#include "moments.hpp"
#include <gtest/gtest.h>
#include <cmath>

using namespace globe::math::polynomial;

TEST(MomentsTest, UnitSphereAreaIsFourPi) {
    Moments moments = Moments::unit_sphere(0);

    EXPECT_NEAR(moments.at(0, 0, 0), 4.0 * M_PI, 1e-12);
}

TEST(MomentsTest, UnitSphereOddMomentsVanish) {
    Moments moments = Moments::unit_sphere(3);

    EXPECT_DOUBLE_EQ(moments.at(1, 0, 0), 0.0);
    EXPECT_DOUBLE_EQ(moments.at(1, 1, 0), 0.0);
    EXPECT_DOUBLE_EQ(moments.at(2, 1, 0), 0.0);
    EXPECT_DOUBLE_EQ(moments.at(1, 1, 1), 0.0);
}

TEST(MomentsTest, UnitSphereSecondMomentsAreFourPiOverThree) {
    Moments moments = Moments::unit_sphere(2);

    EXPECT_NEAR(moments.at(2, 0, 0), 4.0 * M_PI / 3.0, 1e-12);
    EXPECT_NEAR(moments.at(0, 2, 0), 4.0 * M_PI / 3.0, 1e-12);
    EXPECT_NEAR(moments.at(0, 0, 2), 4.0 * M_PI / 3.0, 1e-12);
}

TEST(MomentsTest, UnitSphereFourthMomentsMatchClosedForm) {
    Moments moments = Moments::unit_sphere(4);

    EXPECT_NEAR(moments.at(4, 0, 0), 4.0 * M_PI / 5.0, 1e-12);
    EXPECT_NEAR(moments.at(2, 2, 0), 4.0 * M_PI / 15.0, 1e-12);
}

TEST(MomentsTest, UnitSphereSatisfiesSphereConstraintIdentity) {
    Moments moments = Moments::unit_sphere(6);

    for (const MultiIndex& index : MultiIndex::all_up_to(4)) {
        double sum = 0.0;

        for (int axis = 0; axis < 3; ++axis) {
            sum += moments.at(index.raised(axis).raised(axis));
        }

        EXPECT_NEAR(sum, moments.at(index), 1e-12);
    }
}
