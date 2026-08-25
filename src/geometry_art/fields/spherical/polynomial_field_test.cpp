#include "polynomial_field.hpp"
#include "../../testing/macros.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <vector>

using namespace geometry_art;
using fields::spherical::PolynomialField;
using fields::RegionIntegrals;

namespace {

Polygon octant(int sign_x, int sign_y, int sign_z) {
    VectorS2 x(sign_x, 0, 0);
    VectorS2 y(0, sign_y, 0);
    VectorS2 z(0, 0, sign_z);
    bool counterclockwise = sign_x * sign_y * sign_z > 0;

    if (counterclockwise) {
        return Polygon(std::vector<Arc>{Arc(x, y), Arc(y, z), Arc(z, x)});
    }

    return Polygon(std::vector<Arc>{Arc(x, z), Arc(z, y), Arc(y, x)});
}

std::vector<Polygon> all_octants() {
    std::vector<Polygon> result;

    for (int sign_x : {1, -1}) {
        for (int sign_y : {1, -1}) {
            for (int sign_z : {1, -1}) {
                result.push_back(octant(sign_x, sign_y, sign_z));
            }
        }
    }

    return result;
}

PolynomialField equator_dense_field() {
    Eigen::Matrix3d quadratic_form = Eigen::Matrix3d::Zero();
    quadratic_form(2, 2) = -0.9;
    return PolynomialField::quadratic(1.0, Vector3(0.1, -0.2, 0.3), quadratic_form);
}

}

TEST(PolynomialFieldTest, SecondMomentOfConstantFieldOverEquatorialQuarterArc) {
    PolynomialField field = PolynomialField::constant(1.0);
    Arc quarter(VectorS2(1, 0, 0), VectorS2(0, 1, 0));

    Matrix3 moment = field.second_moment(quarter);

    EXPECT_NEAR(moment(0, 0), M_PI / 4.0, 1e-12);
    EXPECT_NEAR(moment(1, 1), M_PI / 4.0, 1e-12);
    EXPECT_NEAR(moment(0, 1), 0.5, 1e-12);
    EXPECT_NEAR(moment(1, 0), 0.5, 1e-12);
    EXPECT_NEAR(moment(2, 2), 0.0, 1e-12);
    EXPECT_NEAR(moment(0, 2), 0.0, 1e-12);
    EXPECT_NEAR(moment(1, 2), 0.0, 1e-12);
}

TEST(PolynomialFieldTest, SecondMomentTraceIsTheArcMass) {
    PolynomialField field = PolynomialField::linear(2.0, Vector3(0.3, -0.2, 0.5));
    Arc arc(VectorS2(1, 0, 0), VectorS2(0, 1, 1).normalized());

    Matrix3 moment = field.second_moment(arc);

    EXPECT_NEAR(moment.trace(), field.integrals(arc).mass, 1e-12);
}

// For a linear density the gradient is constant, so each partial's second
// moment is that component times the arc's bare second moment -- a closed
// form the arc already computes.
TEST(PolynomialFieldTest, GradientSecondMomentsOfLinearFieldAreClosedForm) {
    Vector3 gradient(0.4, -0.3, 0.8);
    PolynomialField field = PolynomialField::linear(1.0, gradient);
    Arc arc(VectorS2(1, 0, 0), VectorS2(0, 1, 1).normalized());

    std::array<Matrix3, 3> moments = field.gradient_second_moments(arc);
    Matrix3 bare = arc.second_moment();

    for (int axis = 0; axis < 3; ++axis) {
        EXPECT_LT((moments[axis] - gradient[axis] * bare).norm(), 1e-12);
    }
}

TEST(PolynomialFieldTest, GradientSecondMomentsOfConstantFieldVanish) {
    PolynomialField field = PolynomialField::constant(2.0);
    Arc arc(VectorS2(1, 0, 0), VectorS2(0, 0, 1));

    for (const Matrix3& moment : field.gradient_second_moments(arc)) {
        EXPECT_EQ(moment.norm(), 0.0);
    }
}

TEST(PolynomialFieldTest, ConstantFieldValueAndTotalMass) {
    PolynomialField field = PolynomialField::constant(2.0);

    EXPECT_DOUBLE_EQ(field.value(VectorS2(0, 0, 1)), 2.0);
    EXPECT_NEAR(field.total_mass(), 8.0 * M_PI, 1e-12);
}

TEST(PolynomialFieldTest, LinearFieldValue) {
    PolynomialField field = PolynomialField::linear(1.0, Vector3(0, 0, 2.0));

    EXPECT_DOUBLE_EQ(field.value(VectorS2(0, 0, 1)), 3.0);
    EXPECT_DOUBLE_EQ(field.value(VectorS2(1, 0, 0)), 1.0);
}

TEST(PolynomialFieldTest, QuadraticTotalMassIncludesTraceTerm) {
    PolynomialField field = equator_dense_field();

    EXPECT_NEAR(field.total_mass(), 4.0 * M_PI * (1.0 - 0.3), 1e-12);
}

TEST(PolynomialFieldTest, OctantMassOfLinearFieldIsClosedForm) {
    PolynomialField field = PolynomialField::linear(1.0, Vector3(0, 0, 2.0));

    RegionIntegrals integrals = field.integrals(octant(1, 1, 1));

    EXPECT_NEAR(integrals.mass, M_PI / 2.0 + 2.0 * M_PI / 4.0, 1e-12);
}

TEST(PolynomialFieldTest, OctantFirstMomentOfLinearFieldIsClosedForm) {
    PolynomialField field = PolynomialField::linear(1.0, Vector3(0, 0, 2.0));

    RegionIntegrals integrals = field.integrals(octant(1, 1, 1));

    EXPECT_NEAR(integrals.first_moment.x(), M_PI / 4.0 + 2.0 / 3.0, 1e-12);
    EXPECT_NEAR(integrals.first_moment.y(), M_PI / 4.0 + 2.0 / 3.0, 1e-12);
    EXPECT_NEAR(integrals.first_moment.z(), M_PI / 4.0 + 2.0 * M_PI / 6.0, 1e-12);
}

TEST(PolynomialFieldTest, QuarterArcIntegralsOfLinearFieldAreClosedForm) {
    PolynomialField field = PolynomialField::linear(1.0, Vector3(2.0, 3.0, 0.0));

    RegionIntegrals integrals = field.integrals(Arc(VectorS2(1, 0, 0), VectorS2(0, 1, 0)));

    EXPECT_NEAR(integrals.mass, M_PI / 2.0 + 2.0 + 3.0, 1e-12);
    EXPECT_NEAR(integrals.first_moment.x(), 1.0 + 2.0 * M_PI / 4.0 + 3.0 * 0.5, 1e-12);
    EXPECT_NEAR(integrals.first_moment.y(), 1.0 + 2.0 * 0.5 + 3.0 * M_PI / 4.0, 1e-12);
    EXPECT_NEAR(integrals.first_moment.z(), 0.0, 1e-12);
}

TEST(PolynomialFieldTest, OctantsPartitionTotalMassAndFirstMoment) {
    PolynomialField field = equator_dense_field();
    double mass = 0.0;
    Vector3 first_moment = Vector3::Zero();

    for (const Polygon& piece : all_octants()) {
        RegionIntegrals integrals = field.integrals(piece);
        mass += integrals.mass;
        first_moment += integrals.first_moment;
    }

    Vector3 expected_first_moment = Vector3(0.1, -0.2, 0.3) * 4.0 * M_PI / 3.0;
    EXPECT_NEAR(mass, field.total_mass(), 1e-12);
    EXPECT_NEAR((first_moment - expected_first_moment).norm(), 0.0, 1e-12);
}

TEST(PolynomialFieldTest, EXPENSIVE_OctantIntegralsMatchQuadratureOfValue) {
    REQUIRE_EXPENSIVE();

    PolynomialField field = equator_dense_field();
    Polygon region = octant(1, 1, 1);
    RegionIntegrals integrals = field.integrals(region);

    constexpr size_t STEPS = 2000;
    double mass = 0.0;
    Vector3 first_moment = Vector3::Zero();

    for (size_t i = 0; i < STEPS; ++i) {
        double z = (i + 0.5) / STEPS;
        double ring = std::sqrt(1.0 - z * z);

        for (size_t j = 0; j < STEPS; ++j) {
            double phi = (M_PI / 2.0) * (j + 0.5) / STEPS;
            VectorS2 point(ring * std::cos(phi), ring * std::sin(phi), z);
            double weight = (M_PI / 2.0) / (STEPS * STEPS);
            mass += field.value(point) * weight;
            first_moment += field.value(point) * weight * point;
        }
    }

    EXPECT_NEAR(integrals.mass, mass, 1e-6);
    EXPECT_NEAR((integrals.first_moment - first_moment).norm(), 0.0, 1e-6);
}
