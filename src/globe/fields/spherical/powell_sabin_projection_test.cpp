#include "powell_sabin_projection.hpp"
#include "piecewise_polynomial_field.hpp"
#include "polynomial_field.hpp"
#include "../scalar/noise_field.hpp"
#include "../../geometry/spherical/powell_sabin_refinement.hpp"
#include "../../geometry/spherical/triangle_mesh.hpp"
#include "../../generators/spherical/fibonacci_point_generator.hpp"
#include "../../math/interval.hpp"
#include "../../testing/macros.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <vector>

using namespace globe;
using fields::spherical::PiecewisePolynomialField;
using fields::spherical::PolynomialField;
using fields::spherical::PowellSabinProjection;
using fields::scalar::NoiseField;
using geometry::spherical::PowellSabinRefinement;
using geometry::spherical::TriangleMesh;

namespace {

struct ConstantScalarField {
    double constant;
    double value(const VectorS2&) const { return constant; }
};

struct SmoothAnalyticField {
    double value(const VectorS2& point) const {
        return 0.6 + 0.3 * std::sin(3.0 * point.x()) * std::cos(2.0 * point.y()) + 0.1 * point.z() * point.z();
    }
};

struct QuadraticScalarField {
    PolynomialField reference;
    double value(const VectorS2& point) const { return reference.value(point); }
};

// No linear term: a piece is homogeneous of degree two, and on the sphere a
// linear function is a harmonic of a different degree, so no such piece can
// hold one.
PolynomialField tilted_quadratic() {
    Eigen::Matrix3d form;
    form << 0.7, 0.2, -0.1,
            0.2, -0.4, 0.3,
           -0.1, 0.3, 0.5;

    return PolynomialField::quadratic(1.0, Vector3::Zero(), form);
}

// The derivative along a great circle, read from one side of a point only,
// so that a kink shows up as a disagreement between the two sides rather
// than being averaged away.
template<typename FieldType>
double one_sided_derivative(const FieldType& field, const VectorS2& point, const Vector3& direction, double step) {
    VectorS2 moved = VectorS2(point + step * direction).normalized();
    return (field.value(moved) - field.value(point)) / step;
}

template<typename FieldType>
double worst_kink(const FieldType& field, const TriangleMesh& mesh, double step) {
    double worst = 0.0;

    for (const auto& triangle : mesh.triangles) {
        for (int corner = 0; corner < 3; ++corner) {
            VectorS2 from = mesh.vertices[triangle[corner]];
            VectorS2 to = mesh.vertices[triangle[(corner + 1) % 3]];

            // A point strictly inside the edge, and the direction crossing it.
            VectorS2 seam = VectorS2(0.63 * from + 0.37 * to).normalized();
            Vector3 along = VectorS2(to - from).normalized();
            Vector3 across = Vector3(seam.cross(along)).normalized();

            double ahead = one_sided_derivative(field, seam, across, step);
            double behind = one_sided_derivative(field, seam, Vector3(-across), step);

            worst = std::max(worst, std::abs(ahead + behind));
        }
    }

    return worst;
}

} // namespace

// A piece is a quadratic that is homogeneous in space, so the space it lives
// in holds the constants and the quadratics but not the linear functions: on
// the sphere a linear function is a different harmonic degree entirely.
TEST(PowellSabinProjectionTest, ReproducesAConstantExactly) {
    ConstantScalarField constant{0.7};
    PowellSabinProjection projection;
    auto result = projection.project(TriangleMesh::icosphere(1), constant);

    for (const VectorS2& point : generators::spherical::FibonacciPointGenerator().generate(500)) {
        EXPECT_NEAR(result.field.value(point), 0.7, 1e-12);
    }

    // The bound is honest but not tight: the coefficients of a constant are
    // the constant times the corner dot products, which fall below one.
    EXPECT_GT(result.lowest_coefficient, 0.0);
    EXPECT_LE(result.lowest_coefficient, 0.7 + 1e-12);
}

// A global quadratic is itself a C1 quadratic on the split, so the
// projection has a zero-residual candidate and must return it exactly.
// This is the strongest available check that the smoothness conditions and
// the corner conditions are the right ones.
TEST(PowellSabinProjectionTest, ReproducesAGlobalQuadraticFieldExactly) {
    QuadraticScalarField quadratic{tilted_quadratic()};
    PowellSabinProjection projection;
    auto result = projection.project(TriangleMesh::icosphere(1), quadratic);

    for (const VectorS2& point : generators::spherical::FibonacciPointGenerator().generate(500)) {
        EXPECT_NEAR(result.field.value(point), quadratic.value(point), 1e-9);
    }
}

TEST(PowellSabinProjectionTest, SplitPointsLieInsideTheEdgesTheySplit) {
    PowellSabinRefinement refinement(TriangleMesh::icosphere(2));

    for (const auto& cell : refinement.cells()) {
        for (int split = PowellSabinRefinement::SPLIT_01; split <= PowellSabinRefinement::SPLIT_20; ++split) {
            int from = split - PowellSabinRefinement::SPLIT_01;
            int to = (from + 1) % 3;

            EXPECT_GT(cell.points[split].dot(cell.points[from]), 0.0);
            EXPECT_GT(cell.points[split].dot(cell.points[to]), 0.0);
        }
    }
}

TEST(PowellSabinProjectionTest, RefinementGivesSixTrianglesPerTriangle) {
    TriangleMesh mesh = TriangleMesh::icosphere(2);
    PowellSabinRefinement refinement(mesh);

    EXPECT_EQ(refinement.mesh().triangles.size(), 6 * mesh.triangles.size());
    EXPECT_EQ(refinement.cells().size(), mesh.triangles.size());
}

// The point of the whole construction: the gradient does not jump where the
// pieces meet. The Lagrange field on the same mesh is the control -- it is
// continuous but kinked, and the kink is what makes the constraint curvature
// discontinuous for the optimizer.
TEST(PowellSabinProjectionTest, EXPENSIVE_HasNoKinksWhereTheLagrangeFieldDoes) {
    REQUIRE_EXPENSIVE();
    NoiseField noise(Interval(0.2, 1.0));
    TriangleMesh mesh = TriangleMesh::icosphere(2);

    PowellSabinProjection projection;
    auto smooth = projection.project(mesh, noise);

    NoiseField same_noise(Interval(0.2, 1.0));
    PiecewisePolynomialField kinked = PiecewisePolynomialField::sample(mesh, 2, same_noise);

    constexpr double STEP = 1e-6;
    double smooth_kink = worst_kink(smooth.field, smooth.field.mesh(), STEP);
    double kinked_kink = worst_kink(kinked, kinked.mesh(), STEP);

    EXPECT_LT(smooth_kink, 1e-3);
    EXPECT_GT(kinked_kink, 1e-2);
    EXPECT_LT(smooth_kink, kinked_kink / 100.0);
}

// Positivity is a property of the coefficients, so it is certified rather
// than sampled: no search over the sphere can miss a dip the coefficients
// already rule out. It holds at every resolution because a projection that
// dips below zero has its gradients damped instead of shipped.
TEST(PowellSabinProjectionTest, EXPENSIVE_CertifiesPositivityAtEveryResolution) {
    REQUIRE_EXPENSIVE();

    for (int subdivisions : {1, 2, 3, 4}) {
        NoiseField noise(Interval(0.2, 1.0));
        PowellSabinProjection projection;
        auto result = projection.project(TriangleMesh::icosphere(subdivisions), noise);

        EXPECT_GT(result.lowest_coefficient, 0.0) << "subdivisions " << subdivisions;

        for (const VectorS2& point : generators::spherical::FibonacciPointGenerator().generate(100000)) {
            ASSERT_GE(result.field.value(point), result.lowest_coefficient - 1e-12) << "subdivisions " << subdivisions;
        }
    }
}

// Damping is a fallback, not the normal path: a field the mesh resolves is
// projected as solved, and the report says so.
TEST(PowellSabinProjectionTest, LeavesGradientsAloneWhenTheMeshCarriesThem) {
    QuadraticScalarField quadratic{tilted_quadratic()};
    PowellSabinProjection projection;
    auto result = projection.project(TriangleMesh::icosphere(2), quadratic);

    EXPECT_EQ(result.least_damping, 1.0);
}

// The projection keeps the represented field inside the range of the
// field it came from once the mesh resolves it: what it cannot follow it
// averages, so it does not overshoot below the floor and needs no damping
// to stay positive. That matters beyond positivity -- a density with a
// lower minimum than it was given poses a harder capacity problem, so a
// representation that widens the range changes the question being asked.
TEST(PowellSabinProjectionTest, EXPENSIVE_KeepsTheFloorItWasGiven) {
    REQUIRE_EXPENSIVE();
    NoiseField noise(Interval(0.2, 1.0));
    PowellSabinProjection projection;
    auto result = projection.project(TriangleMesh::icosphere(4), noise);

    EXPECT_GT(result.lowest_coefficient, 0.19);
    EXPECT_EQ(result.least_damping, 1.0);
}

TEST(PowellSabinProjectionTest, EXPENSIVE_ApproachesTheSampledFieldWithRefinement) {
    REQUIRE_EXPENSIVE();
    NoiseField noise(Interval(0.2, 1.0));
    PowellSabinProjection projection;

    NoiseField coarse_noise(Interval(0.2, 1.0));
    auto coarse = projection.project(TriangleMesh::icosphere(2), coarse_noise);
    auto fine = projection.project(TriangleMesh::icosphere(4), noise);

    NoiseField reference(Interval(0.2, 1.0));
    double coarse_error = 0.0;
    double fine_error = 0.0;

    for (const VectorS2& point : generators::spherical::FibonacciPointGenerator().generate(20000)) {
        double expected = reference.value(point);
        coarse_error = std::max(coarse_error, std::abs(coarse.field.value(point) - expected));
        fine_error = std::max(fine_error, std::abs(fine.field.value(point) - expected));
    }

    EXPECT_LT(fine_error, coarse_error / 2.0);
}


// Smoothness is not bought with accuracy: on a target that is itself
// smooth the C1 field tracks it as closely as the Lagrange field on the
// same mesh, and at the same order.
TEST(PowellSabinProjectionTest, EXPENSIVE_MatchesTheLagrangeFieldsOrderOnASmoothTarget) {
    REQUIRE_EXPENSIVE();

    for (int subdivisions : {2, 3, 4}) {
        TriangleMesh mesh = TriangleMesh::icosphere(subdivisions);

        SmoothAnalyticField target;
        PowellSabinProjection projection;
        auto smooth = projection.project(mesh, target);

        SmoothAnalyticField sampled;
        PiecewisePolynomialField lagrange = PiecewisePolynomialField::sample(mesh, 2, sampled);

        double smooth_error = 0.0;
        double lagrange_error = 0.0;

        for (const VectorS2& point : generators::spherical::FibonacciPointGenerator().generate(50000)) {
            double expected = target.value(point);
            smooth_error = std::max(smooth_error, std::abs(smooth.field.value(point) - expected));
            lagrange_error = std::max(lagrange_error, std::abs(lagrange.value(point) - expected));
        }

        EXPECT_LT(smooth_error, 2.0 * lagrange_error) << "subdivisions " << subdivisions;
        EXPECT_EQ(smooth.least_damping, 1.0) << "subdivisions " << subdivisions;
    }
}

// On rough noise the projection beats the Lagrange field on the same
// mesh: interpolating point samples commits to whatever the roughness does
// at the sample points, while the best approximation in the space lets
// what it cannot follow land in the residual. Smoothness costs nothing
// here -- the C1 field is the more accurate one as well.
TEST(PowellSabinProjectionTest, EXPENSIVE_BeatsTheLagrangeFieldOnRoughNoise) {
    REQUIRE_EXPENSIVE();
    TriangleMesh mesh = TriangleMesh::icosphere(4);

    NoiseField smooth_noise(Interval(0.2, 1.0));
    PowellSabinProjection projection;
    auto smooth = projection.project(mesh, smooth_noise);

    NoiseField lagrange_noise(Interval(0.2, 1.0));
    PiecewisePolynomialField lagrange = PiecewisePolynomialField::sample(mesh, 2, lagrange_noise);

    NoiseField reference(Interval(0.2, 1.0));
    double smooth_error = 0.0;
    double lagrange_error = 0.0;

    for (const VectorS2& point : generators::spherical::FibonacciPointGenerator().generate(50000)) {
        double expected = reference.value(point);
        smooth_error = std::max(smooth_error, std::abs(smooth.field.value(point) - expected));
        lagrange_error = std::max(lagrange_error, std::abs(lagrange.value(point) - expected));
    }

    EXPECT_LT(smooth_error, lagrange_error);
}




// What the tessellation actually reads is the density's average over a
// cell, so that is the error worth bounding, and it is the error a
// projection is built to bound: structure finer than the mesh lands in the
// residual, whose local averages are small, rather than being aliased into
// the field. A mesh matched to the cell scale -- the coarsest the
// bandwidth rule permits -- is enough for better than one percent.
TEST(PowellSabinProjectionTest, EXPENSIVE_AveragesWellAtTheCellScale) {
    REQUIRE_EXPENSIVE();
    constexpr int SITES = 200;
    const double cap_radius = 2.0 / std::sqrt(static_cast<double>(SITES));

    auto directions = generators::spherical::FibonacciPointGenerator().generate(400000);
    auto centers = generators::spherical::FibonacciPointGenerator().generate(200);

    auto worst_cell_scale_error = [&](int subdivisions) {
        TriangleMesh mesh = TriangleMesh::icosphere(subdivisions);
        NoiseField target(Interval(0.2, 1.0));
        PowellSabinProjection projection;
        auto result = projection.project(mesh, target);
        NoiseField reference(Interval(0.2, 1.0));
        double worst = 0.0;

        for (const VectorS2& center : centers) {
            double expected = 0.0;
            double actual = 0.0;
            int count = 0;

            for (const VectorS2& direction : directions) {
                if (direction.dot(center) < std::cos(cap_radius)) {
                    continue;
                }

                expected += reference.value(direction);
                actual += result.field.value(direction);
                ++count;
            }

            if (count >= 40) {
                worst = std::max(worst, std::abs(actual - expected) / expected);
            }
        }

        return worst;
    };

    EXPECT_LT(worst_cell_scale_error(3), 0.01);
    EXPECT_LT(worst_cell_scale_error(4), 0.002);
}

