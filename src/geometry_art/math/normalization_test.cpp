#include "normalization.hpp"
#include <Eigen/Geometry>
#include <gtest/gtest.h>
#include <cmath>
#include <vector>

using namespace geometry_art;

namespace {

// A quadratic test function stands in for the energy: its own gradient and
// Hessian are exact, so any disagreement is in the normalisation map.
Matrix3 quadratic_form() {
    Matrix3 form;
    form << 0.7, -0.3, 0.2,
           -0.3, 1.1, 0.4,
            0.2, 0.4, -0.5;
    return form;
}

Vector3 linear_term() {
    return Vector3(0.4, -0.9, 0.25);
}

double site_value(const Vector3& site) {
    return 0.5 * site.dot(quadratic_form() * site) + linear_term().dot(site);
}

double normalized_value(const Vector3& point) {
    return site_value(point.normalized());
}

Vector3 site_gradient(const Vector3& site) {
    return quadratic_form() * site + linear_term();
}

double second_difference(const Vector3& point, int row, int column, double step) {
    Vector3 row_step = Vector3::Unit(row) * step;
    Vector3 column_step = Vector3::Unit(column) * step;

    return (
        normalized_value(point + row_step + column_step) -
        normalized_value(point + row_step - column_step) -
        normalized_value(point - row_step + column_step) +
        normalized_value(point - row_step - column_step)
    ) / (4.0 * step * step);
}

} // namespace

TEST(NormalizationTest, GradientMatchesFiniteDifferences) {
    Vector3 point(0.6, -0.4, 0.9);
    Normalization normalization(point);
    Vector3 gradient = normalization.gradient(site_gradient(normalization.site()));
    double step = 1e-6;

    for (int axis = 0; axis < 3; ++axis) {
        Vector3 offset = Vector3::Unit(axis) * step;
        double expected = (normalized_value(point + offset) - normalized_value(point - offset)) / (2.0 * step);
        EXPECT_NEAR(gradient(axis), expected, 1e-7);
    }
}

TEST(NormalizationTest, HessianMatchesFiniteDifferences) {
    Vector3 point(0.6, -0.4, 0.9);
    Normalization normalization(point);
    Matrix3 hessian = normalization.hessian(site_gradient(normalization.site()), quadratic_form());
    double step = 1e-5;

    for (int row = 0; row < 3; ++row) {
        for (int column = 0; column < 3; ++column) {
            EXPECT_NEAR(hessian(row, column), second_difference(point, row, column, step), 1e-5);
        }
    }
}

TEST(NormalizationTest, HessianAnnihilatesTheRadialDirectionUpToTheGradient) {
    Vector3 point(0.6, -0.4, 0.9);
    Normalization normalization(point);
    Vector3 gradient = normalization.gradient(site_gradient(normalization.site()));
    Matrix3 hessian = normalization.hessian(site_gradient(normalization.site()), quadratic_form());

    Vector3 radial = hessian * point + gradient;

    EXPECT_NEAR(radial.norm(), 0.0, 1e-12);
}

TEST(NormalizationTest, TangentialHessianHasNoRadialComponent) {
    Vector3 point(0.6, -0.4, 0.9);
    Normalization normalization(point);
    Matrix3 tangential = normalization.tangential_hessian(site_gradient(normalization.site()), quadratic_form());

    EXPECT_NEAR((tangential * normalization.site()).norm(), 0.0, 1e-12);
    EXPECT_NEAR((normalization.site().transpose() * tangential).norm(), 0.0, 1e-12);
}

TEST(NormalizationTest, TangentialHessianAgreesWithTheAmbientOneOnTangentDirections) {
    Vector3 point(0.6, -0.4, 0.9);
    Normalization normalization(point);
    Vector3 gradient = site_gradient(normalization.site());
    Matrix3 ambient = normalization.hessian(gradient, quadratic_form());
    Matrix3 tangential = normalization.tangential_hessian(gradient, quadratic_form());

    Vector3 axis(0.0, 0.0, 1.0);
    Vector3 direction = Vector3(normalization.site().cross(axis)).normalized();

    EXPECT_NEAR(direction.dot(ambient * direction), direction.dot(tangential * direction), 1e-12);
}
