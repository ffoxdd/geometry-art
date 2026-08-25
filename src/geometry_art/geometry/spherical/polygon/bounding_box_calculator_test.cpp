#include "gtest/gtest.h"
#include "bounding_box_calculator.hpp"
#include "../arc.hpp"
#include "../../../testing/assertions/geometric.hpp"
#include <cmath>

using namespace geometry_art;

TEST(PolygonBoundingBoxCalculatorTest, ArcCrossingEquator) {
    std::vector<Arc> arcs{
        Arc(VectorS2(0.5, 0, std::sqrt(0.75)), VectorS2(0, -1, 0), VectorS2(0, 1, 0)),
        Arc(VectorS2(0, -1, 0), VectorS2(-0.5, 0, std::sqrt(0.75)), VectorS2(0, 1, 0)),
        Arc(VectorS2(-0.5, 0, std::sqrt(0.75)), VectorS2(0.5, 0, std::sqrt(0.75)), VectorS2(0, 1, 0)),
    };

    auto bounding_box = PolygonBoundingBoxCalculator(arcs).calculate();

    EXPECT_LE(bounding_box.z_interval().low(), 0.0);
    EXPECT_GE(bounding_box.z_interval().high(), std::sqrt(0.75));
}

TEST(PolygonBoundingBoxCalculatorTest, HorizontalArc) {
    std::vector<Arc> arcs{
        Arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0), VectorS2(0, 0, 1)),
        Arc(VectorS2(0, 1, 0), VectorS2(1, 0, 0), VectorS2(0, 0, 1)),
    };

    auto bounding_box = PolygonBoundingBoxCalculator(arcs).calculate();

    EXPECT_NEAR(bounding_box.z_interval().low(), 0.0, 1e-9);
    EXPECT_NEAR(bounding_box.z_interval().high(), 0.0, 1e-9);
}


TEST(PolygonBoundingBoxCalculatorTest, SmallArcNearEquator) {
    const double eps = 0.1;
    const double z_low = std::sqrt(1 - eps * eps);

    std::vector<Arc> arcs{
        Arc(VectorS2(eps, 0, z_low), VectorS2(0, eps, z_low), VectorS2(0, 0, 1)),
        Arc(VectorS2(0, eps, z_low), VectorS2(eps, 0, z_low), VectorS2(0, 0, 1)),
    };

    auto bounding_box = PolygonBoundingBoxCalculator(arcs).calculate();

    EXPECT_GE(bounding_box.z_interval().low(), z_low - 1e-6);
    EXPECT_LE(bounding_box.z_interval().high(), z_low + 1e-6);
}

TEST(PolygonBoundingBoxCalculatorTest, LargeArcAlmostFullCircle) {
    const double sqrt2 = std::sqrt(2.0);
    const double eps = 0.01;

    std::vector<Arc> arcs{
        Arc(VectorS2(sqrt2/2, sqrt2/2, 0), VectorS2(-sqrt2/2 + eps, sqrt2/2, 0), VectorS2(0, 0, 1)),
        Arc(VectorS2(-sqrt2/2 + eps, sqrt2/2, 0), VectorS2(sqrt2/2, sqrt2/2, 0), VectorS2(0, 0, 1)),
    };

    auto bounding_box = PolygonBoundingBoxCalculator(arcs).calculate();

    EXPECT_NEAR(bounding_box.z_interval().low(), 0.0, 1e-2);
    EXPECT_NEAR(bounding_box.z_interval().high(), 0.0, 1e-2);
}

TEST(PolygonBoundingBoxCalculatorTest, PolygonWithNorthPoleContained) {
    std::vector<Arc> arcs{
        Arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0), VectorS2(0, 0, 1)),
        Arc(VectorS2(0, 1, 0), VectorS2(-1, 0, 0), VectorS2(0, 0, 1)),
        Arc(VectorS2(-1, 0, 0), VectorS2(0, -1, 0), VectorS2(0, 0, 1)),
        Arc(VectorS2(0, -1, 0), VectorS2(1, 0, 0), VectorS2(0, 0, 1)),
    };

    auto bounding_box = PolygonBoundingBoxCalculator(arcs).calculate();

    EXPECT_NEAR(bounding_box.z_interval().low(), 0.0, 1e-9);
    EXPECT_NEAR(bounding_box.z_interval().high(), 1.0, 1e-9);
    EXPECT_TRUE(bounding_box.theta_interval().is_full());
}

TEST(PolygonBoundingBoxCalculatorTest, PolygonWithSouthPoleContained) {
    std::vector<Arc> arcs{
        Arc(VectorS2(1, 0, 0), VectorS2(0, -1, 0), VectorS2(0, 0, -1)),
        Arc(VectorS2(0, -1, 0), VectorS2(-1, 0, 0), VectorS2(0, 0, -1)),
        Arc(VectorS2(-1, 0, 0), VectorS2(0, 1, 0), VectorS2(0, 0, -1)),
        Arc(VectorS2(0, 1, 0), VectorS2(1, 0, 0), VectorS2(0, 0, -1)),
    };

    auto bounding_box = PolygonBoundingBoxCalculator(arcs).calculate();

    EXPECT_NEAR(bounding_box.z_interval().low(), -1.0, 1e-9);
    EXPECT_NEAR(bounding_box.z_interval().high(), 0.0, 1e-9);
    EXPECT_TRUE(bounding_box.theta_interval().is_full());
}


TEST(PolygonBoundingBoxCalculatorTest, ThetaWrappingAroundZero) {
    double theta1 = 5.8;
    double theta2 = 0.4;
    double z_low = 0.0;
    double z_high = 0.3;
    double r_low = std::sqrt(1.0 - z_low * z_low);
    double r_high = std::sqrt(1.0 - z_high * z_high);

    VectorS2 p1(r_low * std::cos(theta1), r_low * std::sin(theta1), z_low);
    VectorS2 p2(r_low * std::cos(theta2), r_low * std::sin(theta2), z_low);
    VectorS2 p3(r_high * std::cos(theta2), r_high * std::sin(theta2), z_high);
    VectorS2 p4(r_high * std::cos(theta1), r_high * std::sin(theta1), z_high);

    VectorS2 n1 = p1.cross(p2).normalized();
    VectorS2 n2 = p2.cross(p3).normalized();
    VectorS2 n3 = p3.cross(p4).normalized();
    VectorS2 n4 = p4.cross(p1).normalized();

    std::vector<Arc> arcs{
        Arc(p1, p2, n1),
        Arc(p2, p3, n2),
        Arc(p3, p4, n3),
        Arc(p4, p1, n4)
    };

    auto bounding_box = PolygonBoundingBoxCalculator(arcs).calculate();

    EXPECT_GT(bounding_box.theta_interval().start(), M_PI);
    EXPECT_GT(bounding_box.theta_interval().end(), TWO_PI);
}
