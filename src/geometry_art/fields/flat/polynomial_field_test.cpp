#include "polynomial_field.hpp"
#include "../../geometry/planar/domain.hpp"
#include "../../geometry/planar/polygon.hpp"
#include "../../geometry/planar/segment.hpp"
#include <gtest/gtest.h>

using namespace geometry_art;
using fields::flat::PolynomialField;
using geometry::planar::Domain;
using geometry::planar::Polygon;
using geometry::planar::Segment;

TEST(FlatPolynomialFieldTest, ConstantIntegratesARectangleExactly) {
    PolynomialField field = PolynomialField::constant(2.0, Domain::torus(3.0, 1.0));
    Polygon rectangle = Polygon::rectangle(Vector2(1.0, 2.0), Vector2(3.0, 5.0));

    auto integrals = field.integrals(rectangle);

    EXPECT_NEAR(integrals.mass, 12.0, 1e-12);
    EXPECT_NEAR(integrals.first_moment.x(), 12.0 * 2.0, 1e-12);
    EXPECT_NEAR(integrals.first_moment.y(), 12.0 * 3.5, 1e-12);
    EXPECT_NEAR(integrals.first_moment.z(), 0.0, 1e-12);
    EXPECT_NEAR(field.total_mass(), 6.0, 1e-12);
}

TEST(FlatPolynomialFieldTest, ConstantIntegratesASegmentByLength) {
    PolynomialField field = PolynomialField::constant(0.5, Domain::torus(1.0, 1.0));
    Segment segment(Vector2(0.0, 0.0), Vector2(3.0, 4.0));

    auto integrals = field.integrals(segment);

    EXPECT_NEAR(integrals.mass, 2.5, 1e-12);
    EXPECT_NEAR(integrals.first_moment.x(), 2.5 * 1.5, 1e-12);
    EXPECT_NEAR(integrals.first_moment.y(), 2.5 * 2.0, 1e-12);
    EXPECT_NEAR(field.gradient_masses(segment).norm(), 0.0, 1e-12);
}

TEST(FlatPolynomialFieldTest, ConstantSquaredNormMomentOnTheUnitSquare) {
    PolynomialField field = PolynomialField::constant(3.0, Domain::plane(1.0, 1.0));
    Polygon square = Polygon::rectangle(Vector2(0.0, 0.0), Vector2(1.0, 1.0));

    // integral of x^2 + y^2 over the unit square is 2/3.
    EXPECT_NEAR(field.squared_norm_moment(square), 2.0, 1e-12);
}

TEST(FlatPolynomialFieldTest, ConstantSecondMomentOfAnAxisSegment) {
    PolynomialField field = PolynomialField::constant(1.0, Domain::plane(1.0, 1.0));
    Segment segment(Vector2(0.0, 1.0), Vector2(2.0, 1.0));

    Matrix3 moment = field.second_moment(segment);

    EXPECT_NEAR(moment(0, 0), 8.0 / 3.0, 1e-12);
    EXPECT_NEAR(moment(0, 1), 2.0, 1e-12);
    EXPECT_NEAR(moment(1, 1), 2.0, 1e-12);
    EXPECT_NEAR(moment(2, 2), 0.0, 1e-12);
}

// rho = 1 + 2x on a 3 x 1 plane: mass 3 + 9 = 12.
TEST(FlatPolynomialFieldTest, LinearTotalMassIsTheClosedFormValue) {
    PolynomialField field = PolynomialField::linear(1.0, Vector2(2.0, 0.0), Domain::plane(3.0, 1.0));

    EXPECT_EQ(field.degree(), 1);
    EXPECT_NEAR(field.total_mass(), 12.0, 1e-12);
    EXPECT_NEAR(field.value(Vector2(3.0, 0.5)), 7.0, 1e-12);
}

TEST(FlatPolynomialFieldTest, LinearIntegratesARectangleExactly) {
    PolynomialField field = PolynomialField::linear(1.0, Vector2(2.0, 0.0), Domain::plane(3.0, 1.0));
    Polygon rectangle = Polygon::rectangle(Vector2(1.0, 0.0), Vector2(2.0, 1.0));

    // mass = integral of 1 + 2x over [1,2] = 1 + 3 = 4; x moment = 1.5 + 2 * 7/3.
    auto integrals = field.integrals(rectangle);

    EXPECT_NEAR(integrals.mass, 4.0, 1e-12);
    EXPECT_NEAR(integrals.first_moment.x(), 1.5 + 14.0 / 3.0, 1e-12);
    EXPECT_NEAR(integrals.first_moment.y(), 2.0, 1e-12);
}

TEST(FlatPolynomialFieldTest, LinearGradientIntegralsFollowTheSlope) {
    PolynomialField field = PolynomialField::linear(1.0, Vector2(2.0, 0.0), Domain::plane(3.0, 1.0));
    Segment segment(Vector2(1.0, 0.0), Vector2(1.0, 0.5));

    Vector3 gradient_mass = field.gradient_masses(segment);
    Matrix3 gradient_first = field.gradient_first_moments(segment);

    EXPECT_NEAR(gradient_mass.x(), 2.0 * 0.5, 1e-12);
    EXPECT_NEAR(gradient_mass.y(), 0.0, 1e-12);
    EXPECT_NEAR(gradient_first(0, 0), 2.0 * 0.5 * 1.0, 1e-12);
    EXPECT_NEAR(gradient_first(1, 0), 2.0 * 0.125, 1e-12);
    EXPECT_NEAR(gradient_first(0, 1), 0.0, 1e-12);
}

// A gradient up the cylinder's walled height is constant around its
// wrapped width, which is all periodicity asks of it.
TEST(FlatPolynomialFieldTest, LinearUpAWalledAxisIsAdmittedBesideAWrappedOne) {
    PolynomialField field = PolynomialField::linear(1.0, Vector2(0.0, 3.0), Domain::cylinder(2.0, 1.0));

    EXPECT_NEAR(field.total_mass(), 2.0 + 3.0, 1e-12);
}
