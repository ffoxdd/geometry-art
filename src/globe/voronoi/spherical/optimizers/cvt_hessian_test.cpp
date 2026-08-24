#include "cvt_hessian.hpp"
#include "capacity_constrained_lagrangian.hpp"
#include "../../../fields/spherical/polynomial_field.hpp"
#include "../../../math/normalization.hpp"
#include "../../../testing/macros.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <memory>
#include <vector>

using namespace globe;
using namespace globe::voronoi;
using namespace globe::voronoi::spherical;
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

std::unique_ptr<Sphere> build_sphere(const std::vector<VectorS2>& sites) {
    auto sphere = std::make_unique<Sphere>();

    for (const VectorS2& site : sites) {
        sphere->insert(cgal::to_point(site.normalized()));
    }

    return sphere;
}

PolynomialField test_field() {
    Matrix3 quadratic_form = Matrix3::Zero();
    quadratic_form(2, 2) = -0.6;
    quadratic_form(0, 1) = 0.2;
    quadratic_form(1, 0) = 0.2;
    return PolynomialField::quadratic(1.0, Vector3(0.3, -0.2, 0.1), quadratic_form);
}

// The Lagrangian with no multipliers and no penalty is the bare CVT energy,
// so its site gradients are exactly the derivatives the Hessian differentiates.
std::vector<Vector3> cvt_gradients(const Sphere& sphere, const PolynomialField& field) {
    CapacityConstrainedLagrangian<PolynomialField> lagrangian(field, field.total_mass() / sphere.size());
    return lagrangian.evaluate(sphere, std::vector<double>(sphere.size(), 0.0), 0.0).site_gradients;
}

// The gradient the optimizer actually descends: the CVT gradient carried
// through the normalisation of each ambient site vector.
std::vector<Vector3> normalized_gradients(
    const std::vector<Vector3>& points,
    const PolynomialField& field
) {
    std::vector<VectorS2> sites;
    sites.reserve(points.size());

    for (const Vector3& point : points) {
        sites.push_back(VectorS2(point.normalized()));
    }

    auto sphere = build_sphere(sites);
    std::vector<Vector3> site_gradients = cvt_gradients(*sphere, field);
    std::vector<Vector3> gradients(points.size());

    for (size_t k = 0; k < points.size(); ++k) {
        gradients[k] = Normalization(points[k]).gradient(site_gradients[k]);
    }

    return gradients;
}

std::vector<Vector3> ambient_points(const std::vector<VectorS2>& sites) {
    std::vector<Vector3> points;
    points.reserve(sites.size());

    for (const VectorS2& site : sites) {
        points.emplace_back(site.x(), site.y(), site.z());
    }

    return points;
}

} // namespace

TEST(CvtHessianTest, RowsSumToZero) {
    auto sphere = build_sphere(fibonacci_sites(12));
    HessianBlocks blocks = CvtHessian<PolynomialField>(test_field()).assemble(*sphere);

    for (size_t k = 0; k < blocks.diagonal.size(); ++k) {
        Matrix3 row = blocks.diagonal[k];

        for (const NeighborBlock& block : blocks.neighbors[k]) {
            row += block.value;
        }

        EXPECT_NEAR(row.norm(), 0.0, 1e-12);
    }
}

TEST(CvtHessianTest, BlocksAreSymmetricAcrossSharedBisectors) {
    auto sphere = build_sphere(fibonacci_sites(12));
    HessianBlocks blocks = CvtHessian<PolynomialField>(test_field()).assemble(*sphere);

    for (size_t k = 0; k < blocks.neighbors.size(); ++k) {
        for (const NeighborBlock& block : blocks.neighbors[k]) {
            const std::vector<NeighborBlock>& opposite = blocks.neighbors[block.neighbor_index];
            auto match = std::find_if(
                opposite.begin(),
                opposite.end(),
                [k](const NeighborBlock& candidate) { return candidate.neighbor_index == k; }
            );

            ASSERT_NE(match, opposite.end());
            EXPECT_NEAR((block.value - match->value.transpose()).norm(), 0.0, 1e-12);
        }
        EXPECT_NEAR((blocks.diagonal[k] - blocks.diagonal[k].transpose()).norm(), 0.0, 1e-12);
    }
}

TEST(CvtHessianTest, EXPENSIVE_MatchesFiniteDifferencesOfTheNormalizedGradient) {
    REQUIRE_EXPENSIVE();

    PolynomialField field = test_field();
    std::vector<VectorS2> sites = fibonacci_sites(12);
    std::vector<Vector3> points = ambient_points(sites);
    auto sphere = build_sphere(sites);

    HessianBlocks blocks = CvtHessian<PolynomialField>(field).assemble(*sphere);
    std::vector<Vector3> site_gradients = cvt_gradients(*sphere, field);
    std::vector<Normalization> normalizations;

    for (const Vector3& point : points) {
        normalizations.emplace_back(point);
    }

    double step = 1e-6;

    for (size_t column = 0; column < points.size(); ++column) {
        for (int axis = 0; axis < 3; ++axis) {
            std::vector<Vector3> forward = points;
            std::vector<Vector3> backward = points;
            forward[column] += Vector3::Unit(axis) * step;
            backward[column] -= Vector3::Unit(axis) * step;

            std::vector<Vector3> high = normalized_gradients(forward, field);
            std::vector<Vector3> low = normalized_gradients(backward, field);

            for (size_t row = 0; row < points.size(); ++row) {
                Vector3 expected = (high[row] - low[row]) / (2.0 * step);
                Matrix3 block = Matrix3::Zero();

                if (row == column) {
                    block = normalizations[row].hessian(site_gradients[row], blocks.diagonal[row]);
                } else {
                    for (const NeighborBlock& neighbor : blocks.neighbors[row]) {
                        if (neighbor.neighbor_index == column) {
                            block = normalizations[row].mixed_hessian(neighbor.value, normalizations[column]);
                        }
                    }
                }

                Vector3 predicted = block * Vector3::Unit(axis);
                EXPECT_NEAR((predicted - expected).norm(), 0.0, 1e-4)
                    << "row " << row << " column " << column << " axis " << axis;
            }
        }
    }
}
