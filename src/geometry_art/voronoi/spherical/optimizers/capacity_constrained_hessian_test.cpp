#include "capacity_constrained_hessian.hpp"
#include "cvt_hessian.hpp"
#include "../../../math/normalization.hpp"
#include "capacity_constrained_lagrangian.hpp"
#include "../../../fields/spherical/polynomial_field.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <memory>
#include <vector>

using namespace geometry_art;
using namespace geometry_art::voronoi;
using namespace geometry_art::voronoi::spherical;
using fields::spherical::PolynomialField;

namespace {

std::vector<Vector3> fibonacci_points(size_t count) {
    std::vector<Vector3> points;
    double golden_angle = M_PI * (std::sqrt(5.0) - 1.0);

    for (size_t i = 0; i < count; ++i) {
        double y = 1.0 - 2.0 * (i + 0.5) / count;
        double radius = std::sqrt(1.0 - y * y);
        double theta = golden_angle * i;
        points.emplace_back(std::cos(theta) * radius, y, std::sin(theta) * radius);
    }

    return points;
}

std::unique_ptr<Sphere> build_sphere(const std::vector<Vector3>& points) {
    auto sphere = std::make_unique<Sphere>();

    for (const Vector3& point : points) {
        sphere->insert(cgal::to_point(VectorS2(point.normalized())));
    }

    return sphere;
}

PolynomialField test_field() {
    Matrix3 quadratic_form = Matrix3::Zero();
    quadratic_form(2, 2) = -0.6;
    return PolynomialField::quadratic(1.0, Vector3(0.3, -0.2, 0.1), quadratic_form);
}

struct Fixture {
    std::vector<Vector3> points;
    std::unique_ptr<Sphere> sphere;
    PolynomialField field;
    DiagramState state;
    std::vector<Vector3> site_gradients;

    static Fixture build(size_t count) {
        Fixture fixture{fibonacci_points(count), nullptr, test_field(), {}, {}};
        fixture.sphere = build_sphere(fixture.points);
        CapacityConstrainedLagrangian<PolynomialField> lagrangian(
            fixture.field,
            fixture.field.total_mass() / fixture.sphere->size()
        );
        fixture.state = lagrangian.sphere_state(*fixture.sphere);
        fixture.site_gradients = lagrangian
            .evaluate(*fixture.sphere, std::vector<double>(count, 0.0), 0.0)
            .site_gradients;
        return fixture;
    }

    [[nodiscard]] CapacityConstrainedHessian hessian(double penalty) const {
        return CapacityConstrainedHessian(
            CvtHessian<PolynomialField>(field).assemble(*sphere).template through_manifold<Normalization>(points, site_gradients),
            CapacityJacobian(state),
            points,
            penalty
        );
    }
};

std::vector<Vector3> wave(const std::vector<Vector3>& points, double frequency) {
    std::vector<Vector3> directions(points.size());

    for (size_t k = 0; k < points.size(); ++k) {
        Vector3 raw(std::sin(frequency * k), std::cos(frequency * k + 0.4), std::sin(0.9 * frequency * k));
        directions[k] = raw - raw.dot(points[k]) * points[k];
    }

    return directions;
}

double pairing(const std::vector<Vector3>& left, const std::vector<Vector3>& right) {
    double sum = 0.0;

    for (size_t k = 0; k < left.size(); ++k) {
        sum += left[k].dot(right[k]);
    }

    return sum;
}

} // namespace

TEST(CapacityConstrainedHessianTest, ReducesToTheEnergyCurvatureWithoutPenalty) {
    Fixture fixture = Fixture::build(14);
    HessianBlocks energy = CvtHessian<PolynomialField>(fixture.field)
        .assemble(*fixture.sphere)
        .template through_manifold<Normalization>(fixture.points, fixture.site_gradients);

    std::vector<Vector3> directions = wave(fixture.points, 1.3);
    std::vector<Vector3> combined = fixture.hessian(0.0).multiply(directions);
    std::vector<Vector3> bare = energy.multiply(directions);

    for (size_t k = 0; k < directions.size(); ++k) {
        EXPECT_NEAR((combined[k] - bare[k]).norm(), 0.0, 1e-12);
    }
}

TEST(CapacityConstrainedHessianTest, IsSymmetric) {
    Fixture fixture = Fixture::build(14);
    CapacityConstrainedHessian hessian = fixture.hessian(37.0);

    std::vector<Vector3> left = wave(fixture.points, 1.3);
    std::vector<Vector3> right = wave(fixture.points, 2.7);

    EXPECT_NEAR(pairing(hessian.multiply(left), right), pairing(left, hessian.multiply(right)), 1e-9);
}

TEST(CapacityConstrainedHessianTest, PenaltyStiffensDirectionsThatChangeCapacities) {
    Fixture fixture = Fixture::build(14);
    std::vector<Vector3> directions = wave(fixture.points, 1.3);

    double relaxed = pairing(fixture.hessian(0.0).multiply(directions), directions);
    double stiff = pairing(fixture.hessian(1000.0).multiply(directions), directions);

    EXPECT_GT(stiff, relaxed);
}
