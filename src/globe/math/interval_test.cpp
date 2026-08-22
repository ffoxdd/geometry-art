#include "interval.hpp"
#include <gtest/gtest.h>

using namespace globe;

TEST(IntervalTest, MidpointReturnsMiddleValue) {
    Interval interval(0.0, 10.0);
    EXPECT_DOUBLE_EQ(interval.midpoint(), 5.0);
}

TEST(IntervalTest, MidpointWorksWithNegativeValues) {
    Interval interval(-10.0, 10.0);
    EXPECT_DOUBLE_EQ(interval.midpoint(), 0.0);
}

TEST(IntervalTest, MidpointWorksWithNegativeRange) {
    Interval interval(-20.0, -10.0);
    EXPECT_DOUBLE_EQ(interval.midpoint(), -15.0);
}

TEST(IntervalTest, MidpointOfZeroWidthInterval) {
    Interval interval(5.0, 5.0);
    EXPECT_DOUBLE_EQ(interval.midpoint(), 5.0);
}

TEST(IntervalTest, ClampValueWithinInterval) {
    Interval interval(0.0, 10.0);
    EXPECT_DOUBLE_EQ(interval.clamp(5.0), 5.0);
}

TEST(IntervalTest, ClampValueBelowInterval) {
    Interval interval(0.0, 10.0);
    EXPECT_DOUBLE_EQ(interval.clamp(-5.0), 0.0);
}

TEST(IntervalTest, ClampValueAboveInterval) {
    Interval interval(0.0, 10.0);
    EXPECT_DOUBLE_EQ(interval.clamp(15.0), 10.0);
}

TEST(IntervalTest, ClampAtBoundary) {
    Interval interval(0.0, 10.0);
    EXPECT_DOUBLE_EQ(interval.clamp(0.0), 0.0);
    EXPECT_DOUBLE_EQ(interval.clamp(10.0), 10.0);
}
