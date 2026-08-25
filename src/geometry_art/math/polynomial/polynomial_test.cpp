#include "polynomial.hpp"
#include <gtest/gtest.h>
#include <cmath>

using namespace geometry_art::math::polynomial;
using geometry_art::Vector3;

TEST(PolynomialTest, ConstantEvaluatesToItsValue) {
    Polynomial polynomial = Polynomial::constant(2.5);

    EXPECT_DOUBLE_EQ(polynomial.value(Vector3(0.3, -0.2, 0.9)), 2.5);
}

TEST(PolynomialTest, LinearEvaluatesAffineForm) {
    Polynomial polynomial = Polynomial::linear(1.0, Vector3(2.0, 3.0, 4.0));

    EXPECT_DOUBLE_EQ(polynomial.value(Vector3(1.0, 1.0, 1.0)), 10.0);
    EXPECT_DOUBLE_EQ(polynomial.value(Vector3(0.0, 0.0, 0.5)), 3.0);
}

TEST(PolynomialTest, QuadraticEvaluatesQuadraticForm) {
    Eigen::Matrix3d quadratic_form = Eigen::Matrix3d::Zero();
    quadratic_form(0, 1) = 1.0;
    quadratic_form(1, 0) = 1.0;
    quadratic_form(2, 2) = -0.9;

    Polynomial polynomial = Polynomial::quadratic(1.0, Vector3::Zero(), quadratic_form);

    EXPECT_DOUBLE_EQ(polynomial.value(Vector3(0.0, 0.0, 1.0)), 0.1);
    EXPECT_DOUBLE_EQ(polynomial.value(Vector3(1.0, 1.0, 0.0)), 3.0);
}

TEST(PolynomialTest, TimesCoordinateShiftsExponents) {
    Polynomial polynomial = Polynomial::linear(1.0, Vector3(2.0, 0.0, 0.0));

    Polynomial shifted = polynomial.times_coordinate(2);

    EXPECT_EQ(shifted.max_degree(), 2);
    EXPECT_DOUBLE_EQ(shifted.coefficient(MultiIndex{0, 0, 1}), 1.0);
    EXPECT_DOUBLE_EQ(shifted.coefficient(MultiIndex{1, 0, 1}), 2.0);
    EXPECT_DOUBLE_EQ(shifted.coefficient(MultiIndex{0, 0, 0}), 0.0);
}

TEST(PolynomialTest, PartialDerivativeOfQuadraticIsItsGradient) {
    Eigen::Matrix3d form;
    form << 0.5, 0.2, -0.1,
            0.2, -0.3, 0.4,
           -0.1, 0.4, 0.6;

    Polynomial quadratic = Polynomial::quadratic(1.5, Vector3(0.3, -0.7, 0.2), form);
    Vector3 point(0.4, -0.9, 1.3);
    Vector3 expected = Vector3(0.3, -0.7, 0.2) + 2.0 * form * point;

    for (int axis = 0; axis < 3; ++axis) {
        EXPECT_NEAR(quadratic.partial_derivative(axis).value(point), expected[axis], 1e-12);
    }
}

TEST(PolynomialTest, PartialDerivativeOfConstantIsZero) {
    Polynomial constant = Polynomial::constant(4.2);

    for (int axis = 0; axis < 3; ++axis) {
        EXPECT_EQ(constant.partial_derivative(axis).value(Vector3(1.0, 2.0, 3.0)), 0.0);
    }
}

TEST(PolynomialTest, IntegrateAgainstUnitSphereMoments) {
    Eigen::Matrix3d quadratic_form = Eigen::Matrix3d::Zero();
    quadratic_form(2, 2) = -0.9;
    Polynomial polynomial = Polynomial::quadratic(1.0, Vector3(5.0, 5.0, 5.0), quadratic_form);

    double integral = polynomial.integrate(Moments::unit_sphere(2));

    EXPECT_NEAR(integral, 4.0 * M_PI * (1.0 - 0.3), 1e-12);
}

TEST(PolynomialTest, IntegrateAcceptsHigherDegreeMoments) {
    Polynomial polynomial = Polynomial::constant(1.0);

    EXPECT_NEAR(polynomial.integrate(Moments::unit_sphere(4)), 4.0 * M_PI, 1e-12);
}
