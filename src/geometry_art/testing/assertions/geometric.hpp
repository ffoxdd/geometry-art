#ifndef GEOMETRY_ART_TESTING_ASSERTIONS_GEOMETRIC_HPP_
#define GEOMETRY_ART_TESTING_ASSERTIONS_GEOMETRIC_HPP_

#include "../../cgal/types.hpp"
#include "../../geometry/cartesian/bounding_box.hpp"
#include "../../math/interval.hpp"
#include <cmath>
#include <gtest/gtest.h>

namespace geometry_art::testing::assertions {

constexpr double DEFAULT_TOLERANCE = 1e-9;

inline bool is_on_unit_sphere(const cgal::Point3 &point, double tolerance = DEFAULT_TOLERANCE) {
    double distance = std::sqrt(
        point.x() * point.x() +
        point.y() * point.y() +
        point.z() * point.z()
    );

    return std::abs(distance - 1.0) < tolerance;
}

inline bool is_on_unit_sphere(const VectorS2 &point, double tolerance = DEFAULT_TOLERANCE) {
    return std::abs(point.norm() - 1.0) < tolerance;
}

inline bool points_approximately_equal(
    const cgal::Point3 &a,
    const cgal::Point3 &b,
    double tolerance = DEFAULT_TOLERANCE
) {
    double dx = a.x() - b.x();
    double dy = a.y() - b.y();
    double dz = a.z() - b.z();
    double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
    return distance < tolerance;
}

inline void expect_points_equal(
    const cgal::Point3 &actual,
    const cgal::Point3 &expected,
    double tolerance = DEFAULT_TOLERANCE
) {
    EXPECT_TRUE(points_approximately_equal(actual, expected, tolerance))
        << "Expected point (" << expected.x() << ", " << expected.y() << ", " << expected.z() << ") "
        << "but got (" << actual.x() << ", " << actual.y() << ", " << actual.z() << ") "
        << "(distance: " << std::sqrt(
            (actual.x() - expected.x()) * (actual.x() - expected.x()) +
            (actual.y() - expected.y()) * (actual.y() - expected.y()) +
            (actual.z() - expected.z()) * (actual.z() - expected.z())
        ) << ", tolerance: " << tolerance << ")";
}

inline void expect_point_in_box(const cgal::Point3 &point, const BoundingBox &box) {
    EXPECT_GE(point.x(), box.x_interval().low())
        << "Point x-coordinate " << point.x() << " is below box minimum " << box.x_interval().low();

    EXPECT_LE(point.x(), box.x_interval().high())
        << "Point x-coordinate " << point.x() << " is above box maximum " << box.x_interval().high();

    EXPECT_GE(point.y(), box.y_interval().low())
        << "Point y-coordinate " << point.y() << " is below box minimum " << box.y_interval().low();

    EXPECT_LE(point.y(), box.y_interval().high())
        << "Point y-coordinate " << point.y() << " is above box maximum " << box.y_interval().high();

    EXPECT_GE(point.z(), box.z_interval().low())
        << "Point z-coordinate " << point.z() << " is below box minimum " << box.z_interval().low();

    EXPECT_LE(point.z(), box.z_interval().high())
        << "Point z-coordinate " << point.z() << " is above box maximum " << box.z_interval().high();
}

inline void expect_intervals_equal(
    const Interval &actual,
    const Interval &expected,
    double tolerance = DEFAULT_TOLERANCE
) {
    EXPECT_NEAR(actual.low(), expected.low(), tolerance)
        << "Expected interval low " << expected.low() << " but got " << actual.low();
    EXPECT_NEAR(actual.high(), expected.high(), tolerance)
        << "Expected interval high " << expected.high() << " but got " << actual.high();
}

} // namespace geometry_art::testing::assertions

namespace geometry_art::testing {
using assertions::DEFAULT_TOLERANCE;
using assertions::is_on_unit_sphere;
using assertions::points_approximately_equal;
using assertions::expect_points_equal;
using assertions::expect_point_in_box;
using assertions::expect_intervals_equal;
}

#endif //GEOMETRY_ART_TESTING_ASSERTIONS_GEOMETRIC_HPP_
