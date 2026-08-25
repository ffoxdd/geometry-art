#include "polygon.hpp"
#include "../../math/polynomial/polynomial.hpp"
#include "../../testing/macros.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <vector>

using namespace geometry_art;
using geometry_art::math::polynomial::Moments;
using geometry::planar::Polygon;
using geometry::planar::Segment;
using geometry_art::math::polynomial::Polynomial;

namespace {

Polygon unit_square() {
    return Polygon::rectangle(Vector2(0.0, 0.0), Vector2(1.0, 1.0));
}

Polygon right_triangle() {
    return Polygon(std::vector<Vector2>{Vector2(0.0, 0.0), Vector2(1.0, 0.0), Vector2(0.0, 1.0)});
}

} // namespace

TEST(PlanarPolygonTest, AreaOfTheUnitSquare) {
    EXPECT_NEAR(unit_square().area(), 1.0, 1e-14);
}

TEST(PlanarPolygonTest, AreaOfARightTriangle) {
    EXPECT_NEAR(right_triangle().area(), 0.5, 1e-14);
}

TEST(PlanarPolygonTest, ZerothMomentIsTheArea) {
    EXPECT_NEAR(unit_square().moments(0).at(0, 0, 0), 1.0, 1e-14);
    EXPECT_NEAR(right_triangle().moments(0).at(0, 0, 0), 0.5, 1e-14);
}

TEST(PlanarPolygonTest, MomentsOfTheUnitSquareMatchTheClosedForm) {
    Moments moments = unit_square().moments(3);

    // The square is a product domain, so every moment factorises.
    for (int a = 0; a <= 3; ++a) {
        for (int b = 0; a + b <= 3; ++b) {
            double expected = 1.0 / ((a + 1) * (b + 1));
            EXPECT_NEAR(moments.at(a, b, 0), expected, 1e-13) << "a " << a << " b " << b;
        }
    }
}

TEST(PlanarPolygonTest, MomentsOfARightTriangleMatchTheClosedForm) {
    Moments moments = right_triangle().moments(3);

    // Over the corner simplex, the integral of x^a y^b is
    // a! b! / (a + b + 2)!.
    auto factorial = [](int value) {
        double result = 1.0;
        for (int index = 2; index <= value; ++index) {
            result *= index;
        }
        return result;
    };

    for (int a = 0; a <= 3; ++a) {
        for (int b = 0; a + b <= 3; ++b) {
            double expected = factorial(a) * factorial(b) / factorial(a + b + 2);
            EXPECT_NEAR(moments.at(a, b, 0), expected, 1e-13) << "a " << a << " b " << b;
        }
    }
}

TEST(PlanarPolygonTest, MomentsAreTranslationConsistent) {
    Vector2 shift(0.37, -0.82);
    Polygon shifted = Polygon(std::vector<Vector2>{
        Vector2(0.0, 0.0) + shift,
        Vector2(1.0, 0.0) + shift,
        Vector2(0.0, 1.0) + shift
    });

    Moments base = right_triangle().moments(1);
    Moments moved = shifted.moments(1);

    EXPECT_NEAR(moved.at(0, 0, 0), base.at(0, 0, 0), 1e-13);
    EXPECT_NEAR(moved.at(1, 0, 0), base.at(1, 0, 0) + shift.x() * base.at(0, 0, 0), 1e-13);
    EXPECT_NEAR(moved.at(0, 1, 0), base.at(0, 1, 0) + shift.y() * base.at(0, 0, 0), 1e-13);
}

TEST(PlanarPolygonTest, CentroidOfTheUnitSquare) {
    Vector2 centroid = unit_square().centroid();

    EXPECT_NEAR(centroid.x(), 0.5, 1e-13);
    EXPECT_NEAR(centroid.y(), 0.5, 1e-13);
}

TEST(PlanarPolygonTest, ClippingHalvesTheSquare) {
    std::optional<Polygon> clipped = unit_square().clipped_by(Vector2(-1.0, 0.0), Vector2(0.5, 0.0));

    ASSERT_TRUE(clipped.has_value());
    EXPECT_NEAR(clipped->area(), 0.5, 1e-13);
}

TEST(PlanarPolygonTest, ClippingEverythingAwayReportsNothing) {
    EXPECT_FALSE(unit_square().clipped_by(Vector2(-1.0, 0.0), Vector2(-0.5, 0.0)).has_value());
}

TEST(PlanarPolygonTest, ContainsAgreesWithTheBoundary) {
    Polygon square = unit_square();

    EXPECT_TRUE(square.contains(Vector2(0.5, 0.5)));
    EXPECT_FALSE(square.contains(Vector2(1.5, 0.5)));
    EXPECT_FALSE(square.contains(Vector2(0.5, -0.5)));
}

TEST(PlanarSegmentTest, MomentsOfAUnitEdgeAlongX) {
    Segment edge(Vector2(0.0, 0.0), Vector2(1.0, 0.0));
    Moments moments = edge.moments(2);

    EXPECT_NEAR(moments.at(0, 0, 0), 1.0, 1e-14);
    EXPECT_NEAR(moments.at(1, 0, 0), 0.5, 1e-14);
    EXPECT_NEAR(moments.at(2, 0, 0), 1.0 / 3.0, 1e-14);
    EXPECT_NEAR(moments.at(0, 1, 0), 0.0, 1e-14);
}

TEST(PlanarSegmentTest, MomentsScaleWithLength) {
    Segment edge(Vector2(0.0, 0.0), Vector2(3.0, 4.0));

    EXPECT_NEAR(edge.moments(0).at(0, 0, 0), 5.0, 1e-13);
    EXPECT_NEAR(edge.moments(1).at(1, 0, 0), 5.0 * 1.5, 1e-13);
    EXPECT_NEAR(edge.moments(1).at(0, 1, 0), 5.0 * 2.0, 1e-13);
}

TEST(PlanarSegmentTest, NormalPointsOutwardForACounterClockwiseLoop) {
    Segment bottom(Vector2(0.0, 0.0), Vector2(1.0, 0.0));
    Vector2 normal = bottom.normal();

    EXPECT_NEAR(normal.x(), 0.0, 1e-14);
    EXPECT_NEAR(normal.y(), -1.0, 1e-14);
}

// The polynomial layer is shared with the spherical geometry unchanged: a
// planar region reports its moments in the same table, and integrating a
// density is the same contraction.
TEST(PlanarPolygonTest, PolynomialIntegratesOverAPlanarRegion) {
    Polynomial density = Polynomial::linear(2.0, Vector3(3.0, -1.0, 0.0));

    double integral = density.integrate(unit_square().moments(density.max_degree()));

    // Over the unit square: 2 * 1 + 3 * (1/2) - 1 * (1/2).
    EXPECT_NEAR(integral, 3.0, 1e-13);
}
