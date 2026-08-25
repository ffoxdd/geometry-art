#include <gtest/gtest.h>
#include "arc.hpp"
#include "../../testing/macros.hpp"
#include <cmath>

using namespace geometry_art;
using geometry_art::math::polynomial::MultiIndex;

TEST(ArcTest, LengthOfQuarterArc) {
    Arc arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0), VectorS2(0, 0, 1));

    EXPECT_NEAR(arc.length(), M_PI / 2.0, 1e-10);
}

TEST(ArcTest, LengthOfZeroArc) {
    VectorS2 p(1, 0, 0);
    Arc arc(p, p, VectorS2(0, 0, 1));

    EXPECT_NEAR(arc.length(), 0.0, 1e-10);
}

TEST(ArcTest, LengthDoesNotComputeMoments) {
    Arc arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0), VectorS2(0, 0, 1));

    double length1 = arc.length();
    double length2 = arc.length();

    EXPECT_DOUBLE_EQ(length1, length2);
}

TEST(ArcTest, FirstMomentOfQuarterArc) {
    Arc arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0), VectorS2(0, 0, 1));

    VectorS2 moment = arc.first_moment();

    EXPECT_NEAR(moment.x(), 1.0, 1e-10);
    EXPECT_NEAR(moment.y(), 1.0, 1e-10);
    EXPECT_NEAR(moment.z(), 0.0, 1e-10);
}

TEST(ArcTest, FirstMomentTowardsPole) {
    Arc arc(VectorS2(1, 0, 0), VectorS2(0, 0, 1), VectorS2(0, -1, 0));

    VectorS2 moment = arc.first_moment();

    EXPECT_NEAR(moment.x(), 1.0, 1e-10);
    EXPECT_NEAR(moment.y(), 0.0, 1e-10);
    EXPECT_NEAR(moment.z(), 1.0, 1e-10);
}

TEST(ArcTest, SecondMomentIsSymmetric) {
    Arc arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0), VectorS2(0, 0, 1));

    Eigen::Matrix3d moment = arc.second_moment();

    EXPECT_NEAR(moment(0, 1), moment(1, 0), 1e-10);
    EXPECT_NEAR(moment(0, 2), moment(2, 0), 1e-10);
    EXPECT_NEAR(moment(1, 2), moment(2, 1), 1e-10);
}

TEST(ArcTest, SecondMomentTraceEqualsLength) {
    Arc arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0), VectorS2(0, 0, 1));

    double trace = arc.second_moment().trace();

    EXPECT_NEAR(trace, arc.length(), 1e-10);
}

TEST(ArcTest, ZeroArcHasZeroMoments) {
    VectorS2 p(1, 0, 0);
    Arc arc(p, p, VectorS2(0, 0, 1));

    EXPECT_NEAR(arc.first_moment().norm(), 0.0, 1e-10);
    EXPECT_NEAR(arc.second_moment().norm(), 0.0, 1e-10);
}

TEST(ArcTest, ContainsSourceAndTarget) {
    Arc arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0), VectorS2(0, 0, 1));

    EXPECT_TRUE(arc.contains(arc.source()));
    EXPECT_TRUE(arc.contains(arc.target()));
}

TEST(ArcTest, ContainsMidpoint) {
    Arc arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0), VectorS2(0, 0, 1));

    double inv_sqrt2 = 1.0 / std::sqrt(2.0);
    VectorS2 midpoint(inv_sqrt2, inv_sqrt2, 0);

    EXPECT_TRUE(arc.contains(midpoint));
}

TEST(ArcTest, DoesNotContainPointOffArc) {
    Arc arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0), VectorS2(0, 0, 1));

    EXPECT_FALSE(arc.contains(VectorS2(-1, 0, 0)));
    EXPECT_FALSE(arc.contains(VectorS2(0, 0, 1)));
}

TEST(ArcTest, SubarcHasShorterLength) {
    Arc arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0), VectorS2(0, 0, 1));

    double inv_sqrt2 = 1.0 / std::sqrt(2.0);
    VectorS2 midpoint(inv_sqrt2, inv_sqrt2, 0);

    Arc subarc = arc.subarc(midpoint);

    EXPECT_LT(subarc.length(), arc.length());
    EXPECT_NEAR(subarc.length(), arc.length() / 2.0, 1e-10);
}

TEST(ArcTest, InterpolateAtZeroReturnsSource) {
    Arc arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0));

    VectorS2 result = arc.interpolate(0.0);

    EXPECT_NEAR(result.x(), 1.0, 1e-10);
    EXPECT_NEAR(result.y(), 0.0, 1e-10);
    EXPECT_NEAR(result.z(), 0.0, 1e-10);
}

TEST(ArcTest, InterpolateAtOneReturnsTarget) {
    Arc arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0));

    VectorS2 result = arc.interpolate(1.0);

    EXPECT_NEAR(result.x(), 0.0, 1e-10);
    EXPECT_NEAR(result.y(), 1.0, 1e-10);
    EXPECT_NEAR(result.z(), 0.0, 1e-10);
}

TEST(ArcTest, InterpolateAtHalfReturnsMidpoint) {
    Arc arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0));

    VectorS2 result = arc.interpolate(0.5);

    double inv_sqrt2 = 1.0 / std::sqrt(2.0);
    EXPECT_NEAR(result.x(), inv_sqrt2, 1e-10);
    EXPECT_NEAR(result.y(), inv_sqrt2, 1e-10);
    EXPECT_NEAR(result.z(), 0.0, 1e-10);
}

TEST(ArcTest, InterpolatedPointLiesOnArc) {
    Arc arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0));

    EXPECT_TRUE(arc.contains(arc.interpolate(0.25)));
    EXPECT_TRUE(arc.contains(arc.interpolate(0.5)));
    EXPECT_TRUE(arc.contains(arc.interpolate(0.75)));
}

TEST(ArcTest, DoesNotContainPointBehindSourceOnSameGreatCircle) {
    Arc arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0));

    EXPECT_FALSE(arc.contains(VectorS2(std::sqrt(0.5), -std::sqrt(0.5), 0)));
}

TEST(ArcTest, AntipodalArcWithExplicitNormalHasLengthPi) {
    Arc arc(VectorS2(1, 0, 0), VectorS2(-1, 0, 0), VectorS2(0, 0, 1));

    EXPECT_NEAR(arc.length(), M_PI, 1e-12);
}

TEST(ArcTest, AntipodalArcFirstMomentPointsAlongTravelDirection) {
    Arc arc(VectorS2(1, 0, 0), VectorS2(-1, 0, 0), VectorS2(0, 0, 1));

    VectorS2 moment = arc.first_moment();

    EXPECT_NEAR(moment.x(), 0.0, 1e-12);
    EXPECT_NEAR(moment.y(), 2.0, 1e-12);
    EXPECT_NEAR(moment.z(), 0.0, 1e-12);
}

TEST(ArcTest, AntipodalArcInterpolatesThroughTravelDirection) {
    Arc arc(VectorS2(1, 0, 0), VectorS2(-1, 0, 0), VectorS2(0, 0, 1));

    VectorS2 midpoint = arc.interpolate(0.5);

    EXPECT_NEAR(midpoint.y(), 1.0, 1e-12);
}

TEST(ArcTest, MomentsOfQuarterArcMatchClosedFormUpToDegreeThree) {
    Arc arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0));

    auto moments = arc.moments(3);

    EXPECT_NEAR(moments.at(0, 0, 0), M_PI / 2.0, 1e-12);
    EXPECT_NEAR(moments.at(1, 0, 0), 1.0, 1e-12);
    EXPECT_NEAR(moments.at(2, 0, 0), M_PI / 4.0, 1e-12);
    EXPECT_NEAR(moments.at(1, 1, 0), 0.5, 1e-12);
    EXPECT_NEAR(moments.at(3, 0, 0), 2.0 / 3.0, 1e-12);
    EXPECT_NEAR(moments.at(2, 1, 0), 1.0 / 3.0, 1e-12);
    EXPECT_NEAR(moments.at(0, 0, 1), 0.0, 1e-12);
    EXPECT_NEAR(moments.at(1, 0, 2), 0.0, 1e-12);
}

TEST(ArcTest, MomentsOfGeneralArcSatisfySphereConstraintIdentity) {
    Arc arc(VectorS2(0.6, 0.0, 0.8), VectorS2(0.0, -0.6, 0.8).normalized());

    auto moments = arc.moments(6);

    for (const auto& index : MultiIndex::all_up_to(4)) {
        double sum = 0.0;

        for (int axis = 0; axis < 3; ++axis) {
            sum += moments.at(index.raised(axis).raised(axis));
        }

        EXPECT_NEAR(sum, moments.at(index), 1e-12);
    }
}

TEST(ArcTest, SecondMomentAgreesWithGeneralMoments) {
    Arc arc(VectorS2(0.6, 0.0, 0.8), VectorS2(0.0, -0.6, 0.8).normalized());

    Eigen::Matrix3d second = arc.second_moment();
    auto moments = arc.moments(2);

    EXPECT_NEAR(second(0, 2), moments.at(1, 0, 1), 1e-14);
    EXPECT_NEAR(second(1, 1), moments.at(0, 2, 0), 1e-14);
}

TEST(ArcTest, EXPENSIVE_MomentsMatchMidpointQuadrature) {
    REQUIRE_EXPENSIVE();

    Arc arc(VectorS2(0.6, 0.0, 0.8), VectorS2(0.0, -0.6, 0.8).normalized());
    auto moments = arc.moments(4);

    constexpr size_t SAMPLES = 200000;
    double length = arc.length();

    for (const auto& index : MultiIndex::all_up_to(4)) {
        double sum = 0.0;

        for (size_t i = 0; i < SAMPLES; ++i) {
            VectorS2 point = arc.interpolate((i + 0.5) / SAMPLES);
            sum += std::pow(point.x(), index.x) * std::pow(point.y(), index.y) * std::pow(point.z(), index.z);
        }

        EXPECT_NEAR(moments.at(index), sum * length / SAMPLES, 1e-8);
    }
}

TEST(ArcTest, ClippedByKeepsArcEntirelyInsideHalfSpace) {
    Arc arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0));

    auto clipped = arc.clipped_by(VectorS2(1, 1, 0).normalized());

    ASSERT_TRUE(clipped.has_value());
    EXPECT_NEAR(clipped->length(), arc.length(), 1e-12);
}

TEST(ArcTest, ClippedByDropsArcEntirelyOutsideHalfSpace) {
    Arc arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0));

    EXPECT_FALSE(arc.clipped_by(VectorS2(-1, -1, 0).normalized()).has_value());
}

TEST(ArcTest, ClippedByCutsAtGreatCircleCrossing) {
    Arc arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0));

    auto clipped = arc.clipped_by(VectorS2(-1, 1, 0).normalized());

    ASSERT_TRUE(clipped.has_value());
    EXPECT_NEAR(clipped->length(), M_PI / 4.0, 1e-12);
    EXPECT_NEAR((clipped->source() - VectorS2(1, 1, 0).normalized()).norm(), 0.0, 1e-12);
    EXPECT_NEAR((clipped->target() - VectorS2(0, 1, 0)).norm(), 0.0, 1e-12);
}

TEST(ArcTest, ClippedPiecesPartitionTheArc) {
    Arc arc(VectorS2(0.6, 0.0, 0.8), VectorS2(0.0, -0.6, 0.8).normalized());
    VectorS2 normal = VectorS2(0.3, 0.4, -0.2).normalized();

    auto kept = arc.clipped_by(normal);
    auto dropped = arc.clipped_by(-normal);

    ASSERT_TRUE(kept.has_value());
    ASSERT_TRUE(dropped.has_value());
    EXPECT_NEAR(kept->length() + dropped->length(), arc.length(), 1e-12);
    EXPECT_NEAR((kept->first_moment() + dropped->first_moment() - arc.first_moment()).norm(), 0.0, 1e-12);
}
