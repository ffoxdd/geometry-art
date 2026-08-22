#include "polynomial_field_fitter.hpp"
#include "polynomial_field.hpp"
#include "../../generators/spherical/fibonacci_point_generator.hpp"
#include "../../testing/macros.hpp"
#include <gtest/gtest.h>
#include <cmath>

using namespace globe;
using fields::spherical::PolynomialField;
using fields::spherical::PolynomialFieldFitter;
using generators::spherical::FibonacciPointGenerator;

namespace {

struct ExponentialField {
    double value(const VectorS2& point) const { return std::exp(point.z()); }
};

}

TEST(PolynomialFieldFitterTest, RecoversAQuadraticExactly) {
    Eigen::Matrix3d quadratic_form = Eigen::Matrix3d::Zero();
    quadratic_form(0, 1) = 0.5;
    quadratic_form(1, 0) = 0.5;
    quadratic_form(2, 2) = -0.9;
    PolynomialField original = PolynomialField::quadratic(1.0, Vector3(0.2, -0.1, 0.3), quadratic_form);

    PolynomialFieldFitter<> fitter(2, 200, FibonacciPointGenerator());
    auto fit = fitter.fit(original);

    EXPECT_NEAR(fit.root_mean_square_residual, 0.0, 1e-10);
    EXPECT_NEAR(fit.field.value(VectorS2(0.6, 0.0, 0.8)), original.value(VectorS2(0.6, 0.0, 0.8)), 1e-10);
    EXPECT_NEAR(fit.field.total_mass(), original.total_mass(), 1e-10);
}

TEST(PolynomialFieldFitterTest, FitDegreeMatchesRequestedDegree) {
    ExponentialField exponential;

    PolynomialFieldFitter<> fitter(4, 500, FibonacciPointGenerator());
    auto fit = fitter.fit(exponential);

    EXPECT_EQ(fit.field.degree(), 4);
}

TEST(PolynomialFieldFitterTest, ResidualShrinksWithDegree) {
    ExponentialField exponential;

    auto low = PolynomialFieldFitter<>(2, 500, FibonacciPointGenerator()).fit(exponential);
    auto high = PolynomialFieldFitter<>(6, 500, FibonacciPointGenerator()).fit(exponential);

    EXPECT_LT(high.root_mean_square_residual, low.root_mean_square_residual / 10.0);
    EXPECT_LT(high.root_mean_square_residual, 1e-4);
}

TEST(PolynomialFieldFitterTest, FittedTotalMassMatchesQuadrature) {
    ExponentialField exponential;

    auto fit = PolynomialFieldFitter<>(8, 2000, FibonacciPointGenerator()).fit(exponential);

    double expected = 2.0 * M_PI * (std::exp(1.0) - std::exp(-1.0));
    EXPECT_NEAR(fit.field.total_mass(), expected, 1e-4);
}
