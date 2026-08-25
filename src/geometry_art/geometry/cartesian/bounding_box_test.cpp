#include "gtest/gtest.h"
#include "bounding_box.hpp"
#include "../../testing/assertions/geometric.hpp"

using namespace geometry_art;
using geometry_art::testing::expect_intervals_equal;

TEST(BoundingBoxTest, DefaultConstructorCreatesUnitCube) {
    BoundingBox box;

    expect_intervals_equal(box.x_interval(), Interval(-1.0, 1.0));
    expect_intervals_equal(box.y_interval(), Interval(-1.0, 1.0));
    expect_intervals_equal(box.z_interval(), Interval(-1.0, 1.0));
}

TEST(BoundingBoxTest, ConstructorWithExplicitIntervals) {
    Interval x_interval(0.0, 2.0);
    Interval y_interval(-1.0, 1.0);
    Interval z_interval(0.5, 1.5);
    BoundingBox box(x_interval, y_interval, z_interval);

    expect_intervals_equal(box.x_interval(), x_interval);
    expect_intervals_equal(box.y_interval(), y_interval);
    expect_intervals_equal(box.z_interval(), z_interval);
}

TEST(BoundingBoxTest, UnitCubeStaticMethod) {
    BoundingBox unit = BoundingBox::unit_cube();

    expect_intervals_equal(unit.x_interval(), Interval(-1.0, 1.0));
    expect_intervals_equal(unit.y_interval(), Interval(-1.0, 1.0));
    expect_intervals_equal(unit.z_interval(), Interval(-1.0, 1.0));
}

TEST(BoundingBoxTest, UnitCubeMatchesDefaultConstructor) {
    BoundingBox default_box;
    BoundingBox unit_box = BoundingBox::unit_cube();

    expect_intervals_equal(default_box.x_interval(), unit_box.x_interval());
    expect_intervals_equal(default_box.y_interval(), unit_box.y_interval());
    expect_intervals_equal(default_box.z_interval(), unit_box.z_interval());
}

