#include "gtest/gtest.h"
#include "bounding_box.hpp"
#include "../../testing/assertions/geometric.hpp"

using namespace geometry_art;

TEST(SphericalBoundingBoxTest, FullSphereHasCorrectProperties) {
    SphericalBoundingBox box = SphericalBoundingBox::full_sphere();

    EXPECT_TRUE(box.theta_interval().is_full());
    EXPECT_NEAR(box.z_interval().low(), -1.0, 1e-9);
    EXPECT_NEAR(box.z_interval().high(), 1.0, 1e-9);
    EXPECT_NEAR(box.area(), 4 * M_PI, 1e-9);
}

TEST(SphericalBoundingBoxTest, RespectsExplicitIntervals) {
    ThetaInterval theta_interval(M_PI / 6, M_PI / 2);
    Interval z_interval(0.2, 0.5);

    SphericalBoundingBox box(theta_interval, z_interval);

    EXPECT_NEAR(box.theta_interval().start(), theta_interval.start(), 1e-9);
    EXPECT_NEAR(box.theta_interval().measure(), theta_interval.measure(), 1e-9);
    EXPECT_NEAR(box.z_interval().low(), z_interval.low(), 1e-9);
    EXPECT_NEAR(box.z_interval().high(), z_interval.high(), 1e-9);
    EXPECT_DOUBLE_EQ(box.area(), theta_interval.measure() * z_interval.measure());
}

TEST(SphericalBoundingBoxTest, CenterReturnsPointOnUnitSphere) {
    ThetaInterval theta_interval(0.0, M_PI / 2);
    Interval z_interval(0.0, 0.5);
    SphericalBoundingBox box(theta_interval, z_interval);

    VectorS2 center = box.center();

    EXPECT_NEAR(center.norm(), 1.0, 1e-9);
}

TEST(SphericalBoundingBoxTest, CenterUsesIntervalMidpoints) {
    ThetaInterval theta_interval(0.0, M_PI / 2);
    Interval z_interval(-0.5, 0.5);
    SphericalBoundingBox box(theta_interval, z_interval);

    VectorS2 center = box.center();
    double expected_theta = M_PI / 4;
    double expected_z = 0.0;
    double expected_r = 1.0;

    EXPECT_NEAR(center.z(), expected_z, 1e-9);
    EXPECT_NEAR(center.x(), expected_r * std::cos(expected_theta), 1e-9);
    EXPECT_NEAR(center.y(), expected_r * std::sin(expected_theta), 1e-9);
}

TEST(SphericalBoundingBoxTest, HandlesWraparoundThetaInterval) {
    double theta_start = 5.5;
    double theta_measure = (2.0 * M_PI - 5.5) + 0.8;
    ThetaInterval wrapped_theta(theta_start, theta_measure);
    Interval z_interval(0.0, 0.5);
    SphericalBoundingBox box(wrapped_theta, z_interval);

    EXPECT_NEAR(box.theta_interval().measure(), theta_measure, 1e-9);

    VectorS2 center = box.center();
    double expected_theta_center = theta_start + theta_measure / 2.0;
    if (expected_theta_center >= 2.0 * M_PI) {
        expected_theta_center -= 2.0 * M_PI;
    }
    double expected_z = 0.25;
    double expected_r = std::sqrt(1.0 - expected_z * expected_z);

    EXPECT_NEAR(center.z(), expected_z, 1e-9);
    EXPECT_NEAR(center.x(), expected_r * std::cos(expected_theta_center), 1e-9);
    EXPECT_NEAR(center.y(), expected_r * std::sin(expected_theta_center), 1e-9);
}

TEST(SphericalBoundingBoxTest, BoundingSphereRadiusEnclosesBox) {
    ThetaInterval theta_interval(0.0, M_PI / 2);
    Interval z_interval(0.0, 0.5);
    SphericalBoundingBox box(theta_interval, z_interval);

    double radius = box.bounding_sphere_radius();

    double theta_span = box.theta_interval().measure();
    double z_span = box.z_interval().measure();
    double r_max = std::sqrt(1.0 - z_interval.low() * z_interval.low());
    double chord = 2.0 * r_max * std::sin(theta_span / 2.0);
    double expected_radius = std::sqrt(z_span * z_span + chord * chord);

    EXPECT_NEAR(radius, expected_radius, 1e-9);
    EXPECT_GT(radius, 0.0);
}
