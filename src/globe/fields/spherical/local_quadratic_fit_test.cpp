#include "local_quadratic_fit.hpp"
#include "polynomial_field.hpp"
#include "../scalar/noise_field.hpp"
#include "../../math/interval.hpp"
#include "../../testing/macros.hpp"
#include <gtest/gtest.h>
#include <cmath>

using namespace globe;
using fields::spherical::LocalQuadraticFit;
using fields::spherical::PolynomialField;
using fields::scalar::NoiseField;

namespace {

struct ConstantField {
    double constant;
    double value(const VectorS2&) const { return constant; }
};

struct QuadraticField {
    PolynomialField reference;
    double value(const VectorS2& point) const { return reference.value(point); }
};

PolynomialField tilted_quadratic() {
    Eigen::Matrix3d form;
    form << 0.7, 0.2, -0.1,
            0.2, -0.4, 0.3,
           -0.1, 0.3, 0.5;

    return PolynomialField::quadratic(1.0, Vector3::Zero(), form);
}

} // namespace

TEST(LocalQuadraticFitTest, SamplesStayInsideTheRequestedNeighbourhood) {
    LocalQuadraticFit fit(0.2);
    VectorS2 centre = VectorS2(0.3, -0.5, 0.8).normalized();

    for (const VectorS2& sample : fit.samples_around(centre)) {
        EXPECT_NEAR(sample.norm(), 1.0, 1e-12);
        EXPECT_LE(std::acos(std::min(1.0, sample.dot(centre))), 0.2 + 1e-9);
    }
}

TEST(LocalQuadraticFitTest, ReadsAConstantWithNoGradient) {
    ConstantField field{0.7};
    LocalQuadraticFit fit(0.15);
    auto reading = fit.at(field, VectorS2(0.0, 0.0, 1.0));

    EXPECT_NEAR(reading.value, 0.7, 1e-12);
    EXPECT_NEAR(reading.tangential_gradient.norm(), 0.0, 1e-10);
}

// The fit is onto the space the pieces live in, so anything already in that
// space survives it untouched however wide the neighbourhood.
TEST(LocalQuadraticFitTest, ReadsAQuadraticExactly) {
    QuadraticField field{tilted_quadratic()};
    VectorS2 point = VectorS2(0.4, 0.6, -0.7).normalized();

    Matrix3 form;
    form << 0.7, 0.2, -0.1,
            0.2, -0.4, 0.3,
           -0.1, 0.3, 0.5;
    Vector3 exact = 2.0 * form * point;
    exact -= exact.dot(point) * point;

    for (double radius : {0.05, 0.2, 0.5}) {
        LocalQuadraticFit fit(radius);
        auto reading = fit.at(field, point);

        EXPECT_NEAR(reading.value, field.value(point), 1e-9) << "radius " << radius;
        EXPECT_NEAR((reading.tangential_gradient - exact).norm(), 0.0, 1e-8) << "radius " << radius;
    }
}

// The point of it: on a rough field a difference quotient reports the
// roughness, and a fit over a neighbourhood reports the trend. The two
// disagree, and the fit is the one that is stable as the reading moves.
TEST(LocalQuadraticFitTest, EXPENSIVE_IsSteadierThanADifferenceQuotientOnRoughInput) {
    REQUIRE_EXPENSIVE();
    NoiseField noise(Interval(0.2, 1.0));
    VectorS2 centre = VectorS2(0.2, 0.9, 0.3).normalized();
    Vector3 along = Vector3(centre.cross(Vector3::UnitZ())).normalized();

    LocalQuadraticFit fit(0.12);
    double worst_fit_change = 0.0;
    double worst_quotient_change = 0.0;
    Vector3 previous_fit = Vector3::Zero();
    Vector3 previous_quotient = Vector3::Zero();

    for (int step = 0; step <= 40; ++step) {
        VectorS2 point = VectorS2(centre + (0.002 * step) * along).normalized();
        Vector3 from_fit = fit.at(noise, point).tangential_gradient;

        Vector3 helper = std::abs(point.z()) < 0.9 ? Vector3::UnitZ() : Vector3::UnitX();
        Vector3 first = point.cross(helper).normalized();
        Vector3 second = point.cross(first);
        Vector3 from_quotient = Vector3::Zero();

        for (const Vector3& direction : {first, second}) {
            double ahead = noise.value(VectorS2(point + 1e-5 * direction).normalized());
            double behind = noise.value(VectorS2(point - 1e-5 * direction).normalized());
            from_quotient += (ahead - behind) / 2e-5 * direction;
        }

        if (step > 0) {
            worst_fit_change = std::max(worst_fit_change, (from_fit - previous_fit).norm());
            worst_quotient_change = std::max(worst_quotient_change, (from_quotient - previous_quotient).norm());
        }

        previous_fit = from_fit;
        previous_quotient = from_quotient;
    }

    EXPECT_LT(worst_fit_change, worst_quotient_change);
}
