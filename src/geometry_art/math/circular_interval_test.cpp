#include "circular_interval.hpp"
#include <gtest/gtest.h>
#include <vector>
#include <cmath>

using namespace geometry_art;

constexpr double PI = M_PI;
constexpr double TWO_PI = 2.0 * M_PI;
constexpr double HALF_PI = M_PI / 2.0;

using TwoPiInterval = CircularInterval<TWO_PI>;
constexpr double TEN = 10.0;
using TenInterval = CircularInterval<TEN>;

TEST(CircularIntervalTest, ConstructorNormalizesStart) {
    TwoPiInterval interval(TWO_PI + 0.5, 1.0);
    EXPECT_NEAR(interval.start(), 0.5, 1e-10);
    EXPECT_NEAR(interval.measure(), 1.0, 1e-10);
}

TEST(CircularIntervalTest, ConstructorNormalizesNegativeStart) {
    TwoPiInterval interval(-HALF_PI, 1.0);
    EXPECT_NEAR(interval.start(), 1.5 * PI, 1e-10);
}

TEST(CircularIntervalTest, ConstructorClampsMeasureToPeriod) {
    TwoPiInterval interval(0.0, 3.0 * PI);
    EXPECT_NEAR(interval.measure(), TWO_PI, 1e-10);
}

TEST(CircularIntervalTest, EndReturnsSumOfStartAndMeasure) {
    TwoPiInterval interval(HALF_PI, PI);
    EXPECT_NEAR(interval.end(), 1.5 * PI, 1e-10);
}

TEST(CircularIntervalTest, EndCanExceedPeriod) {
    TwoPiInterval interval(1.5 * PI, PI);
    EXPECT_NEAR(interval.end(), 2.5 * PI, 1e-10);
}

TEST(CircularIntervalTest, FromToCreatesCorrectInterval) {
    auto interval = TwoPiInterval::from_to(HALF_PI, PI);
    EXPECT_NEAR(interval.start(), HALF_PI, 1e-10);
    EXPECT_NEAR(interval.measure(), HALF_PI, 1e-10);
}

TEST(CircularIntervalTest, FromToWrapsAround) {
    auto interval = TwoPiInterval::from_to(1.5 * PI, HALF_PI);
    EXPECT_NEAR(interval.start(), 1.5 * PI, 1e-10);
    EXPECT_NEAR(interval.measure(), PI, 1e-10);
}

TEST(CircularIntervalTest, FromToSamePointGivesFullCircle) {
    auto interval = TwoPiInterval::from_to(PI, PI);
    EXPECT_NEAR(interval.measure(), TWO_PI, 1e-10);
}

TEST(CircularIntervalTest, FullReturnsFullCircle) {
    auto interval = TwoPiInterval::full();
    EXPECT_NEAR(interval.start(), 0.0, 1e-10);
    EXPECT_NEAR(interval.measure(), TWO_PI, 1e-10);
    EXPECT_TRUE(interval.is_full());
}

TEST(CircularIntervalTest, IsFullReturnsTrueForFullCircle) {
    TwoPiInterval interval(0.0, TWO_PI);
    EXPECT_TRUE(interval.is_full());
}

TEST(CircularIntervalTest, IsFullReturnsFalseForPartialCircle) {
    TwoPiInterval interval(0.0, PI);
    EXPECT_FALSE(interval.is_full());
}

TEST(CircularIntervalTest, ContainsPointInsideInterval) {
    TwoPiInterval interval(HALF_PI, PI);
    EXPECT_TRUE(interval.contains(PI));
}

TEST(CircularIntervalTest, ContainsPointAtStart) {
    TwoPiInterval interval(HALF_PI, PI);
    EXPECT_TRUE(interval.contains(HALF_PI));
}

TEST(CircularIntervalTest, ContainsPointAtEnd) {
    TwoPiInterval interval(HALF_PI, PI);
    EXPECT_TRUE(interval.contains(1.5 * PI));
}

TEST(CircularIntervalTest, ContainsPointOutsideInterval) {
    TwoPiInterval interval(HALF_PI, PI);
    EXPECT_FALSE(interval.contains(0.0));
}

TEST(CircularIntervalTest, ContainsWrappingIntervalBeforeWrap) {
    TwoPiInterval interval(1.5 * PI, PI);
    EXPECT_TRUE(interval.contains(1.75 * PI));
}

TEST(CircularIntervalTest, ContainsWrappingIntervalAfterWrap) {
    TwoPiInterval interval(1.5 * PI, PI);
    EXPECT_TRUE(interval.contains(0.25 * PI));
}

TEST(CircularIntervalTest, ContainsWrappingIntervalOutside) {
    TwoPiInterval interval(1.5 * PI, PI);
    EXPECT_FALSE(interval.contains(PI));
}

TEST(CircularIntervalTest, FullCircleContainsEverything) {
    auto interval = TwoPiInterval::full();
    EXPECT_TRUE(interval.contains(0.0));
    EXPECT_TRUE(interval.contains(PI));
    EXPECT_TRUE(interval.contains(1.5 * PI));
}

TEST(CircularIntervalTest, HullOfAdjacentIntervals) {
    TwoPiInterval a(0.0, HALF_PI);
    TwoPiInterval b(HALF_PI, HALF_PI);

    auto result = TwoPiInterval::hull(a, b);

    EXPECT_NEAR(result.start(), 0.0, 1e-10);
    EXPECT_NEAR(result.measure(), PI, 1e-10);
}

TEST(CircularIntervalTest, HullOfOverlappingIntervals) {
    TwoPiInterval a(0.0, PI);
    TwoPiInterval b(HALF_PI, PI);

    auto result = TwoPiInterval::hull(a, b);

    EXPECT_NEAR(result.start(), 0.0, 1e-10);
    EXPECT_NEAR(result.measure(), 1.5 * PI, 1e-10);
}

TEST(CircularIntervalTest, HullOfDisjointIntervalsPicksSmallerGap) {
    TwoPiInterval a(0.0, HALF_PI);
    TwoPiInterval b(PI, HALF_PI);

    auto result = TwoPiInterval::hull(a, b);

    EXPECT_NEAR(result.start(), 0.0, 1e-10);
    EXPECT_NEAR(result.measure(), 1.5 * PI, 1e-10);
}

TEST(CircularIntervalTest, HullOfDisjointIntervalsWrapping) {
    TwoPiInterval a(0.0, HALF_PI);
    TwoPiInterval b(1.5 * PI, HALF_PI);

    auto result = TwoPiInterval::hull(a, b);

    EXPECT_NEAR(result.start(), 1.5 * PI, 1e-10);
    EXPECT_NEAR(result.measure(), PI, 1e-10);
}

TEST(CircularIntervalTest, HullWithFullCircleReturnsFullCircle) {
    TwoPiInterval a(0.0, HALF_PI);
    auto full = TwoPiInterval::full();

    auto result = TwoPiInterval::hull(a, full);

    EXPECT_TRUE(result.is_full());
}

TEST(CircularIntervalTest, HullOfOppositeHalvesGivesFullCircle) {
    TwoPiInterval a(0.0, PI);
    TwoPiInterval b(PI, PI);

    auto result = TwoPiInterval::hull(a, b);

    EXPECT_TRUE(result.is_full());
}

TEST(CircularIntervalTest, HullOfContainedIntervalReturnsOuter) {
    TwoPiInterval outer(0.0, PI);
    TwoPiInterval inner(HALF_PI * 0.5, HALF_PI);

    auto result = TwoPiInterval::hull(outer, inner);

    EXPECT_NEAR(result.start(), 0.0, 1e-10);
    EXPECT_NEAR(result.measure(), PI, 1e-10);
}

TEST(CircularIntervalTest, HullOfThreeAdjacentIntervals) {
    std::vector<TwoPiInterval> intervals = {
        TwoPiInterval(0.0, HALF_PI),
        TwoPiInterval(HALF_PI, HALF_PI),
        TwoPiInterval(PI, HALF_PI),
    };

    auto result = TwoPiInterval::hull(intervals);

    EXPECT_NEAR(result.start(), 0.0, 1e-10);
    EXPECT_NEAR(result.measure(), 1.5 * PI, 1e-10);
}


TEST(CircularIntervalTest, HullOfSingleIntervalReturnsSameInterval) {
    std::vector<TwoPiInterval> single = {
        TwoPiInterval(HALF_PI, PI),
    };

    auto result = TwoPiInterval::hull(single);

    EXPECT_NEAR(result.start(), HALF_PI, 1e-10);
    EXPECT_NEAR(result.measure(), PI, 1e-10);
}

TEST(CircularIntervalTest, HullOfIntervalsSpanningFullCircle) {
    std::vector<TwoPiInterval> intervals = {
        TwoPiInterval(0.0, HALF_PI),
        TwoPiInterval(PI, HALF_PI),
        TwoPiInterval(HALF_PI, HALF_PI),
        TwoPiInterval(1.5 * PI, HALF_PI),
    };

    auto result = TwoPiInterval::hull(intervals);

    EXPECT_TRUE(result.is_full());
}

TEST(CircularIntervalTest, CustomPeriodNormalizesCorrectly) {
    TenInterval interval(15.0, 3.0);
    EXPECT_NEAR(interval.start(), 5.0, 1e-10);
    EXPECT_NEAR(interval.measure(), 3.0, 1e-10);
}

TEST(CircularIntervalTest, CustomPeriodContainsWorks) {
    TenInterval interval(8.0, 4.0);
    EXPECT_TRUE(interval.contains(9.0));
    EXPECT_TRUE(interval.contains(1.0));
    EXPECT_FALSE(interval.contains(5.0));
}

TEST(CircularIntervalTest, CustomPeriodHullWorks) {
    TenInterval a(0.0, 2.0);
    TenInterval b(3.0, 2.0);

    auto result = TenInterval::hull(a, b);

    EXPECT_NEAR(result.start(), 0.0, 1e-10);
    EXPECT_NEAR(result.measure(), 5.0, 1e-10);
}

TEST(CircularIntervalTest, CustomPeriodFullWorks) {
    auto interval = TenInterval::full();
    EXPECT_NEAR(interval.measure(), 10.0, 1e-10);
    EXPECT_TRUE(interval.is_full());
    EXPECT_TRUE(interval.contains(0.0));
    EXPECT_TRUE(interval.contains(5.0));
    EXPECT_TRUE(interval.contains(9.9));
}
