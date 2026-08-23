#include "gtest/gtest.h"
#include "polygon.hpp"
#include "../../../testing/macros.hpp"
#include <Eigen/Geometry>
#include "../arc.hpp"
#include <cmath>

using namespace globe;
using globe::math::polynomial::Moments;
using globe::math::polynomial::MultiIndex;

TEST(PolygonTest, SimplePolygon) {
    Polygon spherical_polygon = Polygon(
        std::vector<Arc>{
            Arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0), VectorS2(0, 0, 1)),
            Arc(VectorS2(0, 1, 0), VectorS2(0, 0, 1), VectorS2(1, 0, 0)),
            Arc(VectorS2(0, 0, 1), VectorS2(1, 0, 0), VectorS2(0, 1, 0)),
        }
    );

    const double sqrt3 = std::sqrt(3);

    EXPECT_TRUE(spherical_polygon.contains(VectorS2(1 / sqrt3, 1 / sqrt3, 1 / sqrt3)));
    EXPECT_FALSE(spherical_polygon.contains(VectorS2(-1, 0, 0)));
    EXPECT_FALSE(spherical_polygon.contains(VectorS2(-1 / sqrt3, -1 / sqrt3, -1 / sqrt3)));
}

TEST(PolygonTest, InsideOutPolygon) {
    Polygon spherical_polygon = Polygon(
        std::vector<Arc>{
            Arc(VectorS2(1, 0, 0), VectorS2(0, 0, 1), VectorS2(0, -1, 0)),
            Arc(VectorS2(0, 0, 1), VectorS2(0, 1, 0), VectorS2(-1, 0, 0)),
            Arc(VectorS2(0, 1, 0), VectorS2(1, 0, 0), VectorS2(0, 0, -1)),
        }
    );

    const double sqrt3 = std::sqrt(3);

    EXPECT_FALSE(spherical_polygon.contains(VectorS2(1 / sqrt3, 1 / sqrt3, 1 / sqrt3)));
    EXPECT_TRUE(spherical_polygon.contains(VectorS2(-1, 0, 0)));
    EXPECT_TRUE(spherical_polygon.contains(VectorS2(-1 / sqrt3, -1 / sqrt3, -1 / sqrt3)));
}

TEST(PolygonTest, Hemisphere) {
    Polygon spherical_polygon = Polygon(
        std::vector<Arc>{
            Arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0), VectorS2(0, 0, 1)),
            Arc(VectorS2(0, 1, 0), VectorS2(-1, 0, 0), VectorS2(0, 0, 1)),
            Arc(VectorS2(-1, 0, 0), VectorS2(0, -1, 0), VectorS2(0, 0, 1)),
            Arc(VectorS2(0, -1, 0), VectorS2(1, 0, 0), VectorS2(0, 0, 1)),
        }
    );


    EXPECT_TRUE(spherical_polygon.contains(VectorS2(0, 0, 1)));
    EXPECT_FALSE(spherical_polygon.contains(VectorS2(0, 0, -1)));
}

TEST(PolygonTest, PathologicalHemisphere) {
    Polygon spherical_polygon = Polygon(
        std::vector<Arc>{
            Arc(VectorS2(1, 0, 0), VectorS2(-1, 0, 0), VectorS2(0, 0, 1)),
            Arc(VectorS2(-1, 0, 0), VectorS2(1, 0, 0), VectorS2(0, 0, 1)),
        }
    );

    EXPECT_TRUE(spherical_polygon.contains(VectorS2(0, 0, 1)));
    EXPECT_FALSE(spherical_polygon.contains(VectorS2(0, 0, -1)));
}

TEST(PolygonTest, PointOnArcCircumcircle) {
    Polygon spherical_polygon = Polygon(
        std::vector<Arc>{
            Arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0), VectorS2(0, 0, 1)),
            Arc(VectorS2(0, 1, 0), VectorS2(0, 0, 1), VectorS2(1, 0, 0)),
            Arc(VectorS2(0, 0, 1), VectorS2(1, 0, 0), VectorS2(0, 1, 0)),
        }
    );

    const double sqrt2 = std::sqrt(2);
    VectorS2 point_on_arc(sqrt2 / 2.0, sqrt2 / 2.0, 0);
    EXPECT_TRUE(spherical_polygon.contains(point_on_arc));

    VectorS2 point_on_circle_outside_arc(-sqrt2 / 2.0, sqrt2 / 2.0, 0);
    EXPECT_FALSE(spherical_polygon.contains(point_on_circle_outside_arc));
}

TEST(PolygonTest, PolygonWithWrappedThetaBoundingBox) {
    double theta1 = 5.8;
    double theta2 = 0.4;
    double z_low = 0.0;
    double z_high = 0.4;
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

    Polygon polygon(std::vector<Arc>{
        Arc(p1, p2, n1),
        Arc(p2, p3, n2),
        Arc(p3, p4, n3),
        Arc(p4, p1, n4)
    });

    SphericalBoundingBox bounding_box = polygon.bounding_box();

    EXPECT_GT(bounding_box.theta_interval().end(), TWO_PI);

    VectorS2 center = bounding_box.center();

    double center_theta = std::atan2(center.y(), center.x());
    if (center_theta < 0.0) center_theta += 2.0 * M_PI;

    EXPECT_TRUE(bounding_box.theta_interval().contains(center_theta));
}

TEST(PolygonTest, BoundingBoxWrappedThetaMeasure) {
    double theta1 = 5.8;
    double theta2 = 0.2;
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

    Polygon polygon(std::vector<Arc>{
        Arc(p1, p2, n1),
        Arc(p2, p3, n2),
        Arc(p3, p4, n3),
        Arc(p4, p1, n4)
    });

    SphericalBoundingBox bounding_box = polygon.bounding_box();

    EXPECT_GT(bounding_box.theta_interval().end(), TWO_PI);

    double expected_theta_measure = (2.0 * M_PI - 5.8) + 0.2;
    EXPECT_NEAR(bounding_box.theta_interval().measure(), expected_theta_measure, 0.05);
}

TEST(PolygonTest, BoundingSphereRadiusWithWrappedTheta) {
    double theta_start = 5.5;
    double theta_measure = (TWO_PI - 5.5) + 0.8;
    ThetaInterval wrapped_theta(theta_start, theta_measure);
    Interval z_interval(0.0, 0.5);
    SphericalBoundingBox bounding_box(wrapped_theta, z_interval);

    double radius = bounding_box.bounding_sphere_radius();

    EXPECT_GT(radius, 0.0);

    double theta_span = bounding_box.theta_interval().measure();
    double z_span = bounding_box.z_interval().measure();
    double r_max = std::sqrt(1.0 - z_interval.low() * z_interval.low());
    double chord = 2.0 * r_max * std::sin(theta_span / 2.0);
    double expected_radius = std::sqrt(z_span * z_span + chord * chord);

    EXPECT_NEAR(radius, expected_radius, 1e-9);
}

TEST(PolygonTest, CentroidReturnsPointOnUnitSphere) {
    Polygon polygon = Polygon(
        std::vector<Arc>{
            Arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0), VectorS2(0, 0, 1)),
            Arc(VectorS2(0, 1, 0), VectorS2(0, 0, 1), VectorS2(1, 0, 0)),
            Arc(VectorS2(0, 0, 1), VectorS2(1, 0, 0), VectorS2(0, 1, 0)),
        }
    );

    VectorS2 centroid = polygon.centroid();

    EXPECT_NEAR(centroid.norm(), 1.0, 1e-9);
}

TEST(PolygonTest, CentroidIsOnUnitSphere) {
    Polygon polygon = Polygon(
        std::vector<Arc>{
            Arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0), VectorS2(0, 0, 1)),
            Arc(VectorS2(0, 1, 0), VectorS2(-1, 0, 0), VectorS2(0, 0, 1)),
            Arc(VectorS2(-1, 0, 0), VectorS2(0, -1, 0), VectorS2(0, 0, 1)),
            Arc(VectorS2(0, -1, 0), VectorS2(1, 0, 0), VectorS2(0, 0, 1)),
        }
    );

    VectorS2 centroid = polygon.centroid();

    EXPECT_NEAR(centroid.norm(), 1.0, 1e-9);
}

TEST(PolygonTest, CentroidForSymmetricPolygon) {
    const double sqrt3 = std::sqrt(3);
    Polygon polygon = Polygon(
        std::vector<Arc>{
            Arc(VectorS2(1, 0, 0), VectorS2(0.5, sqrt3 / 2, 0), VectorS2(0, 0, 1)),
            Arc(VectorS2(0.5, sqrt3 / 2, 0), VectorS2(-0.5, sqrt3 / 2, 0), VectorS2(0, 0, 1)),
            Arc(VectorS2(-0.5, sqrt3 / 2, 0), VectorS2(-1, 0, 0), VectorS2(0, 0, 1)),
            Arc(VectorS2(-1, 0, 0), VectorS2(-0.5, -sqrt3 / 2, 0), VectorS2(0, 0, 1)),
            Arc(VectorS2(-0.5, -sqrt3 / 2, 0), VectorS2(0.5, -sqrt3 / 2, 0), VectorS2(0, 0, 1)),
            Arc(VectorS2(0.5, -sqrt3 / 2, 0), VectorS2(1, 0, 0), VectorS2(0, 0, 1)),
        }
    );

    VectorS2 centroid = polygon.centroid();

    EXPECT_NEAR(centroid.norm(), 1.0, 1e-9);
    EXPECT_NEAR(centroid.z(), 1.0, 1e-9);
}

TEST(PolygonTest, ArcThetaExtrema_FullCircle) {
    Polygon polygon(std::vector<Arc>{
        Arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0), VectorS2(0, 0, 1)),
        Arc(VectorS2(0, 1, 0), VectorS2(-1, 0, 0), VectorS2(0, 0, 1)),
        Arc(VectorS2(-1, 0, 0), VectorS2(0, -1, 0), VectorS2(0, 0, 1)),
        Arc(VectorS2(0, -1, 0), VectorS2(1, 0, 0), VectorS2(0, 0, 1))
    });

    ThetaInterval theta_interval = polygon.bounding_box().theta_interval();

    EXPECT_TRUE(theta_interval.is_full());
}

TEST(PolygonTest, BoundingBox_IncludesNorthPoleZ) {
    Polygon polygon(std::vector<Arc>{
        Arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0), VectorS2(0, 0, 1)),
        Arc(VectorS2(0, 1, 0), VectorS2(-1, 0, 0), VectorS2(0, 0, 1)),
        Arc(VectorS2(-1, 0, 0), VectorS2(0, -1, 0), VectorS2(0, 0, 1)),
        Arc(VectorS2(0, -1, 0), VectorS2(1, 0, 0), VectorS2(0, 0, 1))
    });

    EXPECT_TRUE(polygon.contains(VectorS2(0, 0, 1)));

    SphericalBoundingBox bounding_box = polygon.bounding_box();
    Interval z_interval = bounding_box.z_interval();

    EXPECT_LE(z_interval.low(), 0.0);
    EXPECT_GE(z_interval.high(), 1.0);
    EXPECT_TRUE(z_interval.contains(1.0));
}

TEST(PolygonTest, BoundingBox_IncludesSouthPoleZ) {
    Polygon polygon(std::vector<Arc>{
        Arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0), VectorS2(0, 0, -1)),
        Arc(VectorS2(0, 1, 0), VectorS2(-1, 0, 0), VectorS2(0, 0, -1)),
        Arc(VectorS2(-1, 0, 0), VectorS2(0, -1, 0), VectorS2(0, 0, -1)),
        Arc(VectorS2(0, -1, 0), VectorS2(1, 0, 0), VectorS2(0, 0, -1))
    });

    EXPECT_TRUE(polygon.contains(VectorS2(0, 0, -1)));

    SphericalBoundingBox bounding_box = polygon.bounding_box();
    Interval z_interval = bounding_box.z_interval();

    EXPECT_GE(z_interval.high(), 0.0);
    EXPECT_LE(z_interval.low(), -1.0);
    EXPECT_TRUE(z_interval.contains(-1.0));
}

TEST(PolygonTest, BoundingBox_PoleZEvenWhenArcsDontReach) {
    double z_arc = 0.9;
    double r_arc = std::sqrt(1.0 - z_arc * z_arc);

    Polygon polygon(std::vector<Arc>{
        Arc(VectorS2(r_arc, 0, z_arc), VectorS2(0, r_arc, z_arc), VectorS2(0, 0, 1)),
        Arc(VectorS2(0, r_arc, z_arc), VectorS2(-r_arc, 0, z_arc), VectorS2(0, 0, 1)),
        Arc(VectorS2(-r_arc, 0, z_arc), VectorS2(0, -r_arc, z_arc), VectorS2(0, 0, 1)),
        Arc(VectorS2(0, -r_arc, z_arc), VectorS2(r_arc, 0, z_arc), VectorS2(0, 0, 1))
    });

    EXPECT_TRUE(polygon.contains(VectorS2(0, 0, 1)));

    SphericalBoundingBox bounding_box = polygon.bounding_box();
    Interval z_interval = bounding_box.z_interval();

    EXPECT_LE(z_interval.low(), z_arc);
    EXPECT_GE(z_interval.high(), 1.0);
    EXPECT_TRUE(z_interval.contains(1.0));
}

TEST(PolygonTest, BoundingSphereRadiusContainsAllVertices) {
    Polygon polygon(std::vector<Arc>{
        Arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0), VectorS2(0, 0, 1)),
        Arc(VectorS2(0, 1, 0), VectorS2(0, 0, 1), VectorS2(1, 0, 0)),
        Arc(VectorS2(0, 0, 1), VectorS2(1, 0, 0), VectorS2(0, 1, 0)),
    });

    VectorS2 centroid = polygon.centroid();
    double radius = polygon.bounding_sphere_radius();

    for (const auto& point : polygon.points()) {
        double distance = (point - centroid).norm();
        EXPECT_LE(distance, radius + GEOMETRIC_EPSILON);
    }
}

TEST(PolygonTest, BoundingSphereRadiusIsMinimalForSymmetricPolygon) {
    Polygon polygon(std::vector<Arc>{
        Arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0), VectorS2(0, 0, 1)),
        Arc(VectorS2(0, 1, 0), VectorS2(-1, 0, 0), VectorS2(0, 0, 1)),
        Arc(VectorS2(-1, 0, 0), VectorS2(0, -1, 0), VectorS2(0, 0, 1)),
        Arc(VectorS2(0, -1, 0), VectorS2(1, 0, 0), VectorS2(0, 0, 1)),
    });

    double radius = polygon.bounding_sphere_radius();

    double max_distance = 0.0;
    VectorS2 centroid = polygon.centroid();
    for (const auto& point : polygon.points()) {
        double distance = (point - centroid).norm();
        max_distance = std::max(max_distance, distance);
    }

    EXPECT_NEAR(radius, max_distance, GEOMETRIC_EPSILON);
}

TEST(PolygonTest, AreaOfOctant) {
    Polygon polygon(std::vector<Arc>{
        Arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0), VectorS2(0, 0, 1)),
        Arc(VectorS2(0, 1, 0), VectorS2(0, 0, 1), VectorS2(1, 0, 0)),
        Arc(VectorS2(0, 0, 1), VectorS2(1, 0, 0), VectorS2(0, 1, 0)),
    });

    double expected_area = M_PI / 2.0;
    EXPECT_NEAR(polygon.area(), expected_area, 1e-9);
}

TEST(PolygonTest, AreaOfHemisphere) {
    Polygon polygon(std::vector<Arc>{
        Arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0), VectorS2(0, 0, 1)),
        Arc(VectorS2(0, 1, 0), VectorS2(-1, 0, 0), VectorS2(0, 0, 1)),
        Arc(VectorS2(-1, 0, 0), VectorS2(0, -1, 0), VectorS2(0, 0, 1)),
        Arc(VectorS2(0, -1, 0), VectorS2(1, 0, 0), VectorS2(0, 0, 1)),
    });

    double expected_area = 2.0 * M_PI;
    EXPECT_NEAR(polygon.area(), expected_area, 1e-9);
}

TEST(PolygonTest, MomentsArea) {
    Polygon polygon(std::vector<Arc>{
        Arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0), VectorS2(0, 0, 1)),
        Arc(VectorS2(0, 1, 0), VectorS2(0, 0, 1), VectorS2(1, 0, 0)),
        Arc(VectorS2(0, 0, 1), VectorS2(1, 0, 0), VectorS2(0, 1, 0)),
    });

    double expected_area = M_PI / 2.0;
    EXPECT_NEAR(polygon.area(), expected_area, 1e-6);
}

TEST(PolygonTest, MomentsFirstMomentDirection) {
    Polygon polygon(std::vector<Arc>{
        Arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0), VectorS2(0, 0, 1)),
        Arc(VectorS2(0, 1, 0), VectorS2(0, 0, 1), VectorS2(1, 0, 0)),
        Arc(VectorS2(0, 0, 1), VectorS2(1, 0, 0), VectorS2(0, 1, 0)),
    });

    Eigen::Vector3d first_moment = polygon.first_moment();

    EXPECT_GT(first_moment.x(), 0);
    EXPECT_GT(first_moment.y(), 0);
    EXPECT_GT(first_moment.z(), 0);

    EXPECT_NEAR(first_moment.x(), first_moment.y(), 1e-6);
    EXPECT_NEAR(first_moment.y(), first_moment.z(), 1e-6);
}

TEST(PolygonTest, MomentsSecondMomentSymmetric) {
    Polygon polygon(std::vector<Arc>{
        Arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0), VectorS2(0, 0, 1)),
        Arc(VectorS2(0, 1, 0), VectorS2(0, 0, 1), VectorS2(1, 0, 0)),
        Arc(VectorS2(0, 0, 1), VectorS2(1, 0, 0), VectorS2(0, 1, 0)),
    });

    Eigen::Matrix3d second_moment = polygon.second_moment();

    EXPECT_NEAR(second_moment(0, 1), second_moment(1, 0), 1e-10);
    EXPECT_NEAR(second_moment(0, 2), second_moment(2, 0), 1e-10);
    EXPECT_NEAR(second_moment(1, 2), second_moment(2, 1), 1e-10);
}

TEST(PolygonTest, MomentsSecondMomentTrace) {
    Polygon polygon(std::vector<Arc>{
        Arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0), VectorS2(0, 0, 1)),
        Arc(VectorS2(0, 1, 0), VectorS2(0, 0, 1), VectorS2(1, 0, 0)),
        Arc(VectorS2(0, 0, 1), VectorS2(1, 0, 0), VectorS2(0, 1, 0)),
    });

    EXPECT_NEAR(polygon.second_moment().trace(), polygon.area(), 1e-6);
}

namespace {

Polygon octant() {
    return Polygon(std::vector<Arc>{
        Arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0)),
        Arc(VectorS2(0, 1, 0), VectorS2(0, 0, 1)),
        Arc(VectorS2(0, 0, 1), VectorS2(1, 0, 0)),
    });
}

Polygon inside_out_octant() {
    return Polygon(std::vector<Arc>{
        Arc(VectorS2(1, 0, 0), VectorS2(0, 0, 1)),
        Arc(VectorS2(0, 0, 1), VectorS2(0, 1, 0)),
        Arc(VectorS2(0, 1, 0), VectorS2(1, 0, 0)),
    });
}

Polygon northern_hemisphere() {
    return Polygon(std::vector<Arc>{
        Arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0)),
        Arc(VectorS2(0, 1, 0), VectorS2(-1, 0, 0)),
        Arc(VectorS2(-1, 0, 0), VectorS2(0, -1, 0)),
        Arc(VectorS2(0, -1, 0), VectorS2(1, 0, 0)),
    });
}

Polygon irregular_quadrilateral() {
    VectorS2 a = VectorS2(0.9, 0.1, 0.3).normalized();
    VectorS2 b = VectorS2(0.5, 0.7, 0.2).normalized();
    VectorS2 c = VectorS2(0.2, 0.5, 0.8).normalized();
    VectorS2 d = VectorS2(0.6, -0.2, 0.7).normalized();

    return Polygon(std::vector<Arc>{Arc(a, b), Arc(b, c), Arc(c, d), Arc(d, a)});
}

Polygon rotated(const Polygon& polygon, const Eigen::Matrix3d& rotation) {
    std::vector<Arc> arcs;

    for (const Arc& arc : polygon.arcs()) {
        arcs.emplace_back(rotation * arc.source(), rotation * arc.target());
    }

    return Polygon(arcs);
}

}

TEST(PolygonTest, AreaOfInsideOutOctantIsComplement) {
    EXPECT_NEAR(inside_out_octant().area(), 4.0 * M_PI - M_PI / 2.0, 1e-12);
}

TEST(PolygonTest, FirstMomentOfOctantIsClosedForm) {
    VectorS2 moment = octant().first_moment();

    EXPECT_NEAR(moment.x(), M_PI / 4.0, 1e-12);
    EXPECT_NEAR(moment.y(), M_PI / 4.0, 1e-12);
    EXPECT_NEAR(moment.z(), M_PI / 4.0, 1e-12);
}

TEST(PolygonTest, FirstMomentOfInsideOutOctantIsNegated) {
    VectorS2 moment = inside_out_octant().first_moment();

    EXPECT_NEAR(moment.x(), -M_PI / 4.0, 1e-12);
    EXPECT_NEAR(moment.y(), -M_PI / 4.0, 1e-12);
    EXPECT_NEAR(moment.z(), -M_PI / 4.0, 1e-12);
}

TEST(PolygonTest, FirstMomentOfHemisphereIsPiAlongAxis) {
    VectorS2 moment = northern_hemisphere().first_moment();

    EXPECT_NEAR(moment.x(), 0.0, 1e-12);
    EXPECT_NEAR(moment.y(), 0.0, 1e-12);
    EXPECT_NEAR(moment.z(), M_PI, 1e-12);
}

TEST(PolygonTest, SecondMomentOfOctantIsClosedForm) {
    Eigen::Matrix3d moment = octant().second_moment();

    EXPECT_NEAR(moment(0, 0), M_PI / 6.0, 1e-12);
    EXPECT_NEAR(moment(1, 1), M_PI / 6.0, 1e-12);
    EXPECT_NEAR(moment(2, 2), M_PI / 6.0, 1e-12);
    EXPECT_NEAR(moment(0, 1), 1.0 / 3.0, 1e-12);
    EXPECT_NEAR(moment(1, 2), 1.0 / 3.0, 1e-12);
    EXPECT_NEAR(moment(0, 2), 1.0 / 3.0, 1e-12);
}

TEST(PolygonTest, SecondMomentOfHemisphereIsClosedForm) {
    Eigen::Matrix3d moment = northern_hemisphere().second_moment();

    EXPECT_NEAR(moment(0, 0), 2.0 * M_PI / 3.0, 1e-12);
    EXPECT_NEAR(moment(2, 2), 2.0 * M_PI / 3.0, 1e-12);
    EXPECT_NEAR(moment(0, 1), 0.0, 1e-12);
}

TEST(PolygonTest, ThirdMomentsOfOctantAreClosedForm) {
    auto moments = octant().moments(3);

    EXPECT_NEAR(moments.at(3, 0, 0), M_PI / 8.0, 1e-12);
    EXPECT_NEAR(moments.at(1, 1, 1), 1.0 / 8.0, 1e-12);
    EXPECT_NEAR(moments.at(2, 1, 0), M_PI / 16.0, 1e-12);
}

TEST(PolygonTest, OctantsPartitionTheUnitSphereMoments) {
    constexpr int MAX_DEGREE = 5;
    Moments total(MAX_DEGREE);

    for (int sx : {1, -1}) {
        for (int sy : {1, -1}) {
            for (int sz : {1, -1}) {
                Eigen::Matrix3d reflection = Eigen::Vector3d(sx, sy, sz).asDiagonal();
                bool orientation_preserving = sx * sy * sz > 0;
                Polygon piece = rotated(orientation_preserving ? octant() : inside_out_octant(), reflection);

                if (!orientation_preserving) {
                    piece = Polygon(std::vector<Arc>{
                        Arc(reflection * VectorS2(1, 0, 0), reflection * VectorS2(0, 0, 1)),
                        Arc(reflection * VectorS2(0, 0, 1), reflection * VectorS2(0, 1, 0)),
                        Arc(reflection * VectorS2(0, 1, 0), reflection * VectorS2(1, 0, 0)),
                    });
                }

                auto moments = piece.moments(MAX_DEGREE);
                for (const auto& index : MultiIndex::all_up_to(MAX_DEGREE)) {
                    total.add(index, moments.at(index));
                }
            }
        }
    }

    Moments expected = Moments::unit_sphere(MAX_DEGREE);
    for (const auto& index : MultiIndex::all_up_to(MAX_DEGREE)) {
        EXPECT_NEAR(total.at(index), expected.at(index), 1e-11);
    }
}

TEST(PolygonTest, MomentsSatisfySphereConstraintIdentity) {
    auto moments = irregular_quadrilateral().moments(6);

    for (const auto& index : MultiIndex::all_up_to(4)) {
        double sum = 0.0;

        for (int axis = 0; axis < 3; ++axis) {
            sum += moments.at(index.raised(axis).raised(axis));
        }

        EXPECT_NEAR(sum, moments.at(index), 1e-12);
    }
}

TEST(PolygonTest, MomentsAreRotationCovariant) {
    Eigen::Matrix3d rotation = Eigen::AngleAxisd(0.7, Eigen::Vector3d(0.2, -0.5, 0.8).normalized()).toRotationMatrix();
    Polygon original = irregular_quadrilateral();
    Polygon moved = rotated(original, rotation);

    EXPECT_NEAR(moved.area(), original.area(), 1e-12);

    VectorS2 expected_first = rotation * original.first_moment();
    VectorS2 actual_first = moved.first_moment();
    EXPECT_NEAR((expected_first - actual_first).norm(), 0.0, 1e-12);

    Eigen::Matrix3d expected_second = rotation * original.second_moment() * rotation.transpose();
    Eigen::Matrix3d actual_second = moved.second_moment();
    EXPECT_NEAR((expected_second - actual_second).norm(), 0.0, 1e-12);
}

TEST(PolygonTest, CentroidIsNormalizedFirstMoment) {
    Polygon polygon = irregular_quadrilateral();

    VectorS2 centroid = polygon.centroid();

    EXPECT_NEAR((centroid - polygon.first_moment().normalized()).norm(), 0.0, 1e-12);
}

TEST(PolygonTest, EXPENSIVE_MomentsMatchGridQuadrature) {
    REQUIRE_EXPENSIVE();

    Polygon polygon = irregular_quadrilateral();
    auto moments = polygon.moments(4);

    constexpr size_t LATITUDE_STEPS = 2000;
    constexpr size_t LONGITUDE_STEPS = 4000;
    Moments quadrature(4);

    for (size_t i = 0; i < LATITUDE_STEPS; ++i) {
        double z = -1.0 + 2.0 * (i + 0.5) / LATITUDE_STEPS;
        double ring = std::sqrt(1.0 - z * z);

        for (size_t j = 0; j < LONGITUDE_STEPS; ++j) {
            double phi = TWO_PI * (j + 0.5) / LONGITUDE_STEPS;
            VectorS2 point(ring * std::cos(phi), ring * std::sin(phi), z);

            if (!polygon.contains(point)) {
                continue;
            }

            for (const auto& index : MultiIndex::all_up_to(4)) {
                quadrature.add(index, std::pow(point.x(), index.x) * std::pow(point.y(), index.y) * std::pow(point.z(), index.z));
            }
        }
    }

    double weight = 4.0 * M_PI / (LATITUDE_STEPS * LONGITUDE_STEPS);
    for (const auto& index : MultiIndex::all_up_to(4)) {
        EXPECT_NEAR(moments.at(index), quadrature.at(index) * weight, 2e-3);
    }
}

TEST(PolygonTest, ClippedByReturnsWholePolygonWhenInside) {
    Polygon polygon = octant();

    auto clipped = polygon.clipped_by(VectorS2(1, 1, 1).normalized());

    ASSERT_TRUE(clipped.has_value());
    EXPECT_NEAR(clipped->area(), polygon.area(), 1e-12);
}

TEST(PolygonTest, ClippedByReturnsNothingWhenOutside) {
    Polygon polygon = octant();

    EXPECT_FALSE(polygon.clipped_by(VectorS2(-1, -1, -1).normalized()).has_value());
}

TEST(PolygonTest, ClippedByHalvesTheOctantAlongSymmetryPlane) {
    Polygon polygon = octant();

    auto clipped = polygon.clipped_by(VectorS2(1, -1, 0).normalized());

    ASSERT_TRUE(clipped.has_value());
    EXPECT_NEAR(clipped->area(), polygon.area() / 2.0, 1e-12);
    EXPECT_GT(clipped->first_moment().x(), clipped->first_moment().y());
}

TEST(PolygonTest, ClippedPiecesPartitionMomentsOfTheOriginal) {
    Polygon polygon = irregular_quadrilateral();
    VectorS2 normal = VectorS2(0.2, -0.7, 0.4).normalized();

    auto kept = polygon.clipped_by(normal);
    auto dropped = polygon.clipped_by(-normal);

    ASSERT_TRUE(kept.has_value());
    ASSERT_TRUE(dropped.has_value());

    auto kept_moments = kept->moments(3);
    auto dropped_moments = dropped->moments(3);
    auto moments = polygon.moments(3);

    for (const auto& index : MultiIndex::all_up_to(3)) {
        EXPECT_NEAR(kept_moments.at(index) + dropped_moments.at(index), moments.at(index), 1e-12);
    }
}

TEST(PolygonTest, ClippedArcsKeepTheirGreatCirclesAndTheClosingArcTakesThePlanes) {
    Polygon polygon = octant();
    VectorS2 normal = VectorS2(1, -1, 0).normalized();
    auto clipped = polygon.clipped_by(normal);
    ASSERT_TRUE(clipped.has_value());

    size_t closing_arcs = 0;

    for (const Arc& arc : clipped->arcs()) {
        bool on_plane = (arc.normal() - normal).norm() < 1e-12;
        bool on_original = false;

        for (const Arc& original : polygon.arcs()) {
            on_original = on_original || (arc.normal() - original.normal()).norm() < 1e-12;
        }

        closing_arcs += on_plane;
        EXPECT_TRUE(on_plane || on_original);
    }

    EXPECT_EQ(closing_arcs, 1u);
}

TEST(PolygonTest, ClippingThroughAVertexNeighbourhoodIsContinuous) {
    Polygon polygon = octant();
    VectorS2 vertex(0, 0, 1);
    VectorS2 direction = VectorS2(1, -1, 0).normalized();
    std::vector<double> masses;

    for (double offset : {-1e-7, -1e-9, 0.0, 1e-9, 1e-7}) {
        VectorS2 normal = (direction + offset * vertex).normalized();
        auto clipped = polygon.clipped_by(normal);
        masses.push_back(clipped ? clipped->moments(0).at(0, 0, 0) : 0.0);
    }

    for (size_t i = 1; i < masses.size(); ++i) {
        EXPECT_NEAR(masses[i], masses[i - 1], 1e-6);
    }

    EXPECT_NEAR(masses[2], polygon.moments(0).at(0, 0, 0) / 2.0, 1e-9);
}

TEST(PolygonTest, ClippingIsIdempotent) {
    Polygon polygon = irregular_quadrilateral();
    VectorS2 normal = VectorS2(0.2, -0.7, 0.4).normalized();

    auto once = polygon.clipped_by(normal);
    ASSERT_TRUE(once.has_value());
    auto twice = once->clipped_by(normal);
    ASSERT_TRUE(twice.has_value());

    EXPECT_NEAR(twice->area(), once->area(), 1e-12);
}
