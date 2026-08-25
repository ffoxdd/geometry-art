#include "capacity_jacobian.hpp"
#include "spherical/optimizers/capacity_constrained_lagrangian.hpp"
#include "../fields/spherical/polynomial_field.hpp"
#include "../testing/macros.hpp"
#include <Eigen/Geometry>
#include <gtest/gtest.h>
#include <cmath>
#include <memory>
#include <vector>

using namespace geometry_art;
using namespace geometry_art::voronoi;
using namespace geometry_art::voronoi::spherical;
using fields::spherical::PolynomialField;

namespace {

std::vector<VectorS2> fibonacci_sites(size_t count) {
    std::vector<VectorS2> sites;
    double golden_angle = M_PI * (std::sqrt(5.0) - 1.0);

    for (size_t i = 0; i < count; ++i) {
        double y = 1.0 - 2.0 * (i + 0.5) / count;
        double radius = std::sqrt(1.0 - y * y);
        double theta = golden_angle * i;
        sites.emplace_back(std::cos(theta) * radius, y, std::sin(theta) * radius);
    }

    return sites;
}

std::unique_ptr<Sphere> build_sphere(const std::vector<Vector3>& points) {
    auto sphere = std::make_unique<Sphere>();

    for (const Vector3& point : points) {
        sphere->insert(cgal::to_point(VectorS2(point.normalized())));
    }

    return sphere;
}

std::vector<Vector3> ambient_points(const std::vector<VectorS2>& sites) {
    std::vector<Vector3> points;
    points.reserve(sites.size());

    for (const VectorS2& site : sites) {
        points.emplace_back(site.x(), site.y(), site.z());
    }

    return points;
}

PolynomialField test_field() {
    Matrix3 quadratic_form = Matrix3::Zero();
    quadratic_form(2, 2) = -0.6;
    quadratic_form(0, 1) = 0.2;
    quadratic_form(1, 0) = 0.2;
    return PolynomialField::quadratic(1.0, Vector3(0.3, -0.2, 0.1), quadratic_form);
}

CapacityConstrainedLagrangian<PolynomialField> lagrangian_for(const Sphere& sphere, const PolynomialField& field) {
    return CapacityConstrainedLagrangian<PolynomialField>(field, field.total_mass() / sphere.size());
}

std::vector<double> masses(const Sphere& sphere, const PolynomialField& field) {
    std::vector<double> result;

    for (const CellState& cell : lagrangian_for(sphere, field).cell_states(sphere)) {
        result.push_back(cell.mass);
    }

    return result;
}

CapacityJacobian jacobian_for(const Sphere& sphere, const PolynomialField& field) {
    return CapacityJacobian(lagrangian_for(sphere, field).sphere_state(sphere));
}

Vector3 tangent_at(const Vector3& site, int axis) {
    Vector3 reference = Vector3::Unit(axis);
    Vector3 tangent = reference - reference.dot(site) * site;
    return tangent.normalized();
}

} // namespace

TEST(CapacityJacobianTest, ForwardAndTransposeAreAdjoint) {
    PolynomialField field = test_field();
    std::vector<Vector3> points = ambient_points(fibonacci_sites(14));
    auto sphere = build_sphere(points);
    CapacityJacobian jacobian = jacobian_for(*sphere, field);

    std::vector<Vector3> directions(points.size());
    std::vector<double> values(points.size());

    for (size_t k = 0; k < points.size(); ++k) {
        directions[k] = Vector3(std::sin(1.3 * k), std::cos(2.1 * k), std::sin(0.7 * k + 1.0));
        values[k] = std::cos(1.9 * k) + 0.3;
    }

    std::vector<double> forward = jacobian.apply(directions);
    std::vector<Vector3> transpose = jacobian.transpose_apply(values);

    double forward_pairing = 0.0;
    double transpose_pairing = 0.0;

    for (size_t k = 0; k < points.size(); ++k) {
        forward_pairing += forward[k] * values[k];
        transpose_pairing += transpose[k].dot(directions[k]);
    }

    EXPECT_NEAR(forward_pairing, transpose_pairing, 1e-10);
}

TEST(CapacityJacobianTest, ForwardRatesSumToZero) {
    PolynomialField field = test_field();
    std::vector<Vector3> points = ambient_points(fibonacci_sites(14));
    auto sphere = build_sphere(points);
    CapacityJacobian jacobian = jacobian_for(*sphere, field);

    std::vector<Vector3> directions(points.size());

    for (size_t k = 0; k < points.size(); ++k) {
        directions[k] = tangent_at(points[k], static_cast<int>(k % 3));
    }

    std::vector<double> rates = jacobian.apply(directions);
    double total = 0.0;

    for (double rate : rates) {
        total += rate;
    }

    EXPECT_NEAR(total, 0.0, 1e-10);
}

TEST(CapacityJacobianTest, EXPENSIVE_MatchesFiniteDifferencesOfTheCellMasses) {
    REQUIRE_EXPENSIVE();

    PolynomialField field = test_field();
    std::vector<Vector3> points = ambient_points(fibonacci_sites(14));
    auto sphere = build_sphere(points);
    CapacityJacobian jacobian = jacobian_for(*sphere, field);
    double step = 1e-6;

    for (size_t column = 0; column < points.size(); ++column) {
        for (int axis = 0; axis < 3; ++axis) {
            Vector3 tangent = tangent_at(points[column], axis);

            std::vector<Vector3> directions(points.size(), Vector3::Zero());
            directions[column] = tangent;

            std::vector<Vector3> forward = points;
            std::vector<Vector3> backward = points;
            forward[column] = (points[column] + step * tangent).normalized();
            backward[column] = (points[column] - step * tangent).normalized();

            std::vector<double> high = masses(*build_sphere(forward), field);
            std::vector<double> low = masses(*build_sphere(backward), field);
            std::vector<double> predicted = jacobian.apply(directions);

            for (size_t row = 0; row < points.size(); ++row) {
                double expected = (high[row] - low[row]) / (2.0 * step);
                EXPECT_NEAR(predicted[row], expected, 1e-6)
                    << "row " << row << " column " << column << " axis " << axis;
            }
        }
    }
}
