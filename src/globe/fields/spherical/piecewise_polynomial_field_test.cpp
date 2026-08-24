#include "piecewise_polynomial_field.hpp"
#include "polynomial_field.hpp"
#include "../../geometry/spherical/triangle_mesh.hpp"
#include "../../testing/macros.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <vector>

using namespace globe;
using fields::spherical::PiecewisePolynomialField;
using fields::spherical::PolynomialField;
using fields::RegionIntegrals;
using geometry::spherical::TriangleMesh;

namespace {

struct ConstantScalarField {
    double value(const VectorS2&) const { return 2.0; }
};

struct LinearScalarField {
    Vector3 gradient;
    double value(const VectorS2& point) const { return gradient.dot(point); }
};

struct QuadraticScalarField {
    PolynomialField reference;
    double value(const VectorS2& point) const { return reference.value(point); }
};

struct ExponentialScalarField {
    double value(const VectorS2& point) const { return std::exp(point.z()); }
};

Polygon irregular_quadrilateral() {
    VectorS2 a = VectorS2(0.9, 0.1, 0.3).normalized();
    VectorS2 b = VectorS2(0.5, 0.7, 0.2).normalized();
    VectorS2 c = VectorS2(0.2, 0.5, 0.8).normalized();
    VectorS2 d = VectorS2(0.6, -0.2, 0.7).normalized();

    return Polygon(std::vector<Arc>{Arc(a, b), Arc(b, c), Arc(c, d), Arc(d, a)});
}

Polygon octant() {
    return Polygon(std::vector<Arc>{
        Arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0)),
        Arc(VectorS2(0, 1, 0), VectorS2(0, 0, 1)),
        Arc(VectorS2(0, 0, 1), VectorS2(1, 0, 0)),
    });
}

PolynomialField equator_dense_field() {
    Eigen::Matrix3d quadratic_form = Eigen::Matrix3d::Zero();
    quadratic_form(2, 2) = -0.9;
    quadratic_form(0, 1) = 0.3;
    quadratic_form(1, 0) = 0.3;
    return PolynomialField::quadratic(1.0, Vector3::Zero(), quadratic_form);
}

void expect_integrals_match(const RegionIntegrals& actual, const RegionIntegrals& expected, double tolerance) {
    EXPECT_NEAR(actual.mass, expected.mass, tolerance);
    EXPECT_NEAR((actual.first_moment - expected.first_moment).norm(), 0.0, tolerance);
}

}

TEST(PiecewisePolynomialFieldTest, LagrangeNodesOfDegreeTwoAreCornersAndEdgeMidpoints) {
    std::array<VectorS2, 3> corners{VectorS2(1, 0, 0), VectorS2(0, 1, 0), VectorS2(0, 0, 1)};

    auto nodes = PiecewisePolynomialField::lagrange_nodes(corners, 2);

    ASSERT_EQ(nodes.size(), 6u);
    EXPECT_NEAR((nodes[0] - corners[0]).norm(), 0.0, 1e-12);
    EXPECT_NEAR((nodes[1] - VectorS2(1, 1, 0).normalized()).norm(), 0.0, 1e-12);
    EXPECT_NEAR((nodes[5] - corners[2]).norm(), 0.0, 1e-12);
}

TEST(PiecewisePolynomialFieldTest, DegreeTwoReproducesConstantsExactly) {
    ConstantScalarField constant;
    PiecewisePolynomialField field = PiecewisePolynomialField::sample(TriangleMesh::icosphere(2), 2, constant);

    EXPECT_NEAR(field.value(VectorS2(0.3, -0.4, 0.866).normalized()), 2.0, 1e-10);
    EXPECT_NEAR(field.total_mass(), 8.0 * M_PI, 1e-9);
    expect_integrals_match(field.integrals(irregular_quadrilateral()), PolynomialField::constant(2.0).integrals(irregular_quadrilateral()), 1e-9);
}

TEST(PiecewisePolynomialFieldTest, DegreeOneReproducesAGlobalLinearFieldExactly) {
    LinearScalarField linear{Vector3(0.4, -0.3, 0.8)};
    PiecewisePolynomialField field = PiecewisePolynomialField::sample(TriangleMesh::icosphere(2), 1, linear);
    PolynomialField reference = PolynomialField::linear(0.0, linear.gradient);

    expect_integrals_match(field.integrals(irregular_quadrilateral()), reference.integrals(irregular_quadrilateral()), 1e-10);
    EXPECT_NEAR(field.total_mass(), reference.total_mass(), 1e-10);
}

TEST(PiecewisePolynomialFieldTest, DegreeTwoReproducesAGlobalQuadraticFieldExactly) {
    QuadraticScalarField quadratic{equator_dense_field()};
    PiecewisePolynomialField field = PiecewisePolynomialField::sample(TriangleMesh::icosphere(2), 2, quadratic);

    expect_integrals_match(field.integrals(irregular_quadrilateral()), quadratic.reference.integrals(irregular_quadrilateral()), 1e-9);
    EXPECT_NEAR(field.total_mass(), quadratic.reference.total_mass(), 1e-9);
}

TEST(PiecewisePolynomialFieldTest, ArcIntegralsOfAGlobalQuadraticFieldAreExact) {
    QuadraticScalarField quadratic{equator_dense_field()};
    PiecewisePolynomialField field = PiecewisePolynomialField::sample(TriangleMesh::icosphere(2), 2, quadratic);
    Arc arc(VectorS2(0.6, 0.0, 0.8), VectorS2(0.0, -0.6, 0.8).normalized());

    expect_integrals_match(field.integrals(arc), quadratic.reference.integrals(arc), 1e-9);
}

TEST(PiecewisePolynomialFieldTest, IsContinuousAcrossTriangleEdges) {
    ExponentialScalarField exponential;
    PiecewisePolynomialField field = PiecewisePolynomialField::sample(TriangleMesh::icosphere(1), 2, exponential);
    const TriangleMesh& mesh = field.mesh();
    VectorS2 a = mesh.vertices[mesh.triangles[0][0]];
    VectorS2 b = mesh.vertices[mesh.triangles[0][1]];
    Arc edge(a, b);

    VectorS2 on_edge = edge.interpolate(0.3);
    VectorS2 slightly_left = (on_edge + 1e-7 * edge.normal()).normalized();
    VectorS2 slightly_right = (on_edge - 1e-7 * edge.normal()).normalized();

    EXPECT_NEAR(field.value(slightly_left), field.value(slightly_right), 1e-5);
}

TEST(PiecewisePolynomialFieldTest, OctantsPartitionTotalMass) {
    ExponentialScalarField exponential;
    PiecewisePolynomialField field = PiecewisePolynomialField::sample(TriangleMesh::icosphere(2), 2, exponential);
    double sum = 0.0;

    for (int sign_x : {1, -1}) {
        for (int sign_y : {1, -1}) {
            for (int sign_z : {1, -1}) {
                VectorS2 x(sign_x, 0, 0), y(0, sign_y, 0), z(0, 0, sign_z);
                bool counterclockwise = sign_x * sign_y * sign_z > 0;
                Polygon piece = counterclockwise
                    ? Polygon(std::vector<Arc>{Arc(x, y), Arc(y, z), Arc(z, x)})
                    : Polygon(std::vector<Arc>{Arc(x, z), Arc(z, y), Arc(y, x)});
                sum += field.integrals(piece).mass;
            }
        }
    }

    EXPECT_NEAR(sum, field.total_mass(), 1e-9);
}

TEST(PiecewisePolynomialFieldTest, EXPENSIVE_ConvergesToSmoothFieldWithRefinement) {
    REQUIRE_EXPENSIVE();

    ExponentialScalarField exponential;
    Polygon region = octant();
    double exact = (M_PI / 2.0) * (std::exp(1.0) - 1.0);

    double coarse_error = std::abs(PiecewisePolynomialField::sample(TriangleMesh::icosphere(2), 2, exponential).integrals(region).mass - exact);
    double fine_error = std::abs(PiecewisePolynomialField::sample(TriangleMesh::icosphere(4), 2, exponential).integrals(region).mass - exact);

    EXPECT_LT(fine_error, coarse_error / 8.0);
    EXPECT_LT(fine_error, 1e-6);
}
