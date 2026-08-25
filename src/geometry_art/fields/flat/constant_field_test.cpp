#include "constant_field.hpp"
#include "../../geometry/planar/polygon.hpp"
#include "../../geometry/planar/segment.hpp"
#include <gtest/gtest.h>

using namespace geometry_art;
using fields::flat::ConstantField;
using geometry::planar::Polygon;
using geometry::planar::Segment;

TEST(FlatConstantFieldTest, IntegratesARectangleExactly) {
    ConstantField field(2.0, 3.0, 1.0);
    Polygon rectangle = Polygon::rectangle(Vector2(1.0, 2.0), Vector2(3.0, 5.0));

    auto integrals = field.integrals(rectangle);

    EXPECT_NEAR(integrals.mass, 12.0, 1e-12);
    EXPECT_NEAR(integrals.first_moment.x(), 12.0 * 2.0, 1e-12);
    EXPECT_NEAR(integrals.first_moment.y(), 12.0 * 3.5, 1e-12);
    EXPECT_NEAR(integrals.first_moment.z(), 0.0, 1e-12);
    EXPECT_NEAR(field.total_mass(), 6.0, 1e-12);
}

TEST(FlatConstantFieldTest, IntegratesASegmentByLength) {
    ConstantField field(0.5, 1.0, 1.0);
    Segment segment(Vector2(0.0, 0.0), Vector2(3.0, 4.0));

    auto integrals = field.integrals(segment);

    EXPECT_NEAR(integrals.mass, 2.5, 1e-12);
    EXPECT_NEAR(integrals.first_moment.x(), 2.5 * 1.5, 1e-12);
    EXPECT_NEAR(integrals.first_moment.y(), 2.5 * 2.0, 1e-12);
}

TEST(FlatConstantFieldTest, SquaredNormMomentOnTheUnitSquare) {
    ConstantField field(3.0, 1.0, 1.0);
    Polygon square = Polygon::rectangle(Vector2(0.0, 0.0), Vector2(1.0, 1.0));

    // integral of x^2 + y^2 over the unit square is 2/3.
    EXPECT_NEAR(field.squared_norm_moment(square), 2.0, 1e-12);
}

TEST(FlatConstantFieldTest, SecondMomentOfAnAxisSegment) {
    ConstantField field(1.0, 1.0, 1.0);
    Segment segment(Vector2(0.0, 1.0), Vector2(2.0, 1.0));

    Matrix3 moment = field.second_moment(segment);

    EXPECT_NEAR(moment(0, 0), 8.0 / 3.0, 1e-12);
    EXPECT_NEAR(moment(0, 1), 2.0, 1e-12);
    EXPECT_NEAR(moment(1, 1), 2.0, 1e-12);
    EXPECT_NEAR(moment(2, 2), 0.0, 1e-12);
}
