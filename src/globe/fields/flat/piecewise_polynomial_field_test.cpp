#include "piecewise_polynomial_field.hpp"
#include "constant_field.hpp"
#include "noise_field.hpp"
#include "../../geometry/planar/polygon.hpp"
#include "../../geometry/planar/segment.hpp"
#include "../../math/interval.hpp"
#include "../../testing/macros.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <vector>

using namespace globe;
using fields::flat::ConstantField;
using fields::flat::NoiseField;
using fields::flat::PiecewisePolynomialField;
using geometry::planar::Polygon;
using geometry::planar::Segment;
using globe::math::Interval;

namespace {

// A globally quadratic scalar is not periodic, so sampling wraps the nodes
// on the far edges and corrupts exactly the boundary tiles. Away from those
// tiles every answer has a closed form, and the representation holds
// quadratics exactly.
struct QuadraticScalar {
    double value(const Vector2& point) const {
        return 1.0 + 0.5 * point.x() + 0.25 * point.y() * point.y();
    }
};

} // namespace

TEST(FlatPiecewiseFieldTest, ReproducesAQuadraticExactlyInsideTheTile) {
    QuadraticScalar quadratic;
    auto field = PiecewisePolynomialField::sample(1.0, 1.0, 4, 4, quadratic);

    for (double x : {0.1, 0.33, 0.71}) {
        for (double y : {0.05, 0.4, 0.71}) {
            EXPECT_NEAR(field.value(Vector2(x, y)), quadratic.value(Vector2(x, y)), 1e-10);
        }
    }
}

TEST(FlatPiecewiseFieldTest, IntegratesAQuadraticOverAnInteriorPolygonExactly) {
    QuadraticScalar quadratic;
    auto field = PiecewisePolynomialField::sample(1.0, 1.0, 4, 4, quadratic);
    Polygon region = Polygon::rectangle(Vector2(0.2, 0.3), Vector2(0.7, 0.6));

    // mass = integral of 1 + x/2 + y^2/4 over [0.2,0.7]x[0.3,0.6].
    double mass = 0.15 + 0.5 * 0.225 * 0.3 + 0.25 * 0.5 * (0.216 - 0.027) / 3.0;
    auto integrals = field.integrals(region);

    EXPECT_NEAR(integrals.mass, mass, 1e-12);
}

TEST(FlatPiecewiseFieldTest, TotalMassOfASampledConstantIsExact) {
    ConstantField constant(0.7, 2.0, 1.5);
    auto field = PiecewisePolynomialField::sample(2.0, 1.5, 6, 5, constant);

    EXPECT_NEAR(field.total_mass(), 0.7 * 3.0, 1e-12);
}

TEST(FlatPiecewiseFieldTest, WrapsPeriodicallyForProtrudingRegions) {
    QuadraticScalar quadratic;
    auto field = PiecewisePolynomialField::sample(1.0, 1.0, 8, 8, quadratic);

    // The same physical square, once inside the tile and once shifted a
    // full period: identical mass by periodicity of the representation.
    Polygon inside = Polygon::rectangle(Vector2(0.1, 0.1), Vector2(0.3, 0.3));
    Polygon shifted = Polygon::rectangle(Vector2(1.1, -0.9), Vector2(1.3, -0.7));

    EXPECT_NEAR(field.integrals(inside).mass, field.integrals(shifted).mass, 1e-12);

    // The first moment shifts by exactly the period times the mass.
    Vector3 difference = field.integrals(shifted).first_moment - field.integrals(inside).first_moment;
    double mass = field.integrals(inside).mass;

    EXPECT_NEAR(difference.x(), 1.0 * mass, 1e-12);
    EXPECT_NEAR(difference.y(), -1.0 * mass, 1e-12);
}

TEST(FlatPiecewiseFieldTest, SegmentIntegralsMatchQuadrature) {
    QuadraticScalar quadratic;
    auto field = PiecewisePolynomialField::sample(1.0, 1.0, 4, 4, quadratic);
    Segment segment(Vector2(0.15, 0.2), Vector2(0.7, 0.7));

    double sum = 0.0;
    int samples = 20000;

    for (int k = 0; k < samples; ++k) {
        double fraction = (static_cast<double>(k) + 0.5) / samples;
        sum += quadratic.value(segment.interpolate(fraction));
    }

    double quadrature = sum * segment.length() / samples;

    EXPECT_NEAR(field.integrals(segment).mass, quadrature, 1e-6);
}

TEST(FlatPiecewiseFieldTest, SquaredNormMomentShiftsConsistently) {
    QuadraticScalar quadratic;
    auto field = PiecewisePolynomialField::sample(1.0, 1.0, 8, 8, quadratic);

    Polygon inside = Polygon::rectangle(Vector2(0.2, 0.2), Vector2(0.4, 0.5));
    Polygon shifted = Polygon::rectangle(Vector2(1.2, 0.2), Vector2(1.4, 0.5));

    auto integrals = field.integrals(inside);
    double expected = field.squared_norm_moment(inside) +
        2.0 * integrals.first_moment.x() + integrals.mass;

    EXPECT_NEAR(field.squared_norm_moment(shifted), expected, 1e-12);
}

TEST(FlatPiecewiseFieldTest, NoiseSamplingIsContinuousAcrossTheSeams) {
    NoiseField noise(2.0, 1.0, Interval(0.2, 1.0));
    auto field = PiecewisePolynomialField::sample(2.0, 1.0, 16, 8, noise);

    for (double y : {0.13, 0.5, 0.87}) {
        EXPECT_NEAR(field.value(Vector2(0.0, y)), field.value(Vector2(2.0, y)), 1e-10);
    }

    for (double x : {0.2, 1.1, 1.9}) {
        EXPECT_NEAR(field.value(Vector2(x, 0.0)), field.value(Vector2(x, 1.0)), 1e-10);
    }

    EXPECT_GT(field.lowest_sampled_value(), 0.19);
}
