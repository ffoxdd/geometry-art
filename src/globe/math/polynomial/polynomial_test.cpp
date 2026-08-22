#include "polynomial.hpp"
#include <gtest/gtest.h>
#include <cmath>

using namespace globe::math::polynomial;
using globe::Vector3;

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
