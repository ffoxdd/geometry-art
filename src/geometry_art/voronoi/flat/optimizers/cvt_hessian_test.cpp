#include "cvt_hessian.hpp"
#include "capacity_constrained_lagrangian.hpp"
#include "../core/diagram.hpp"
#include "../../hessian_blocks.hpp"
#include "../../../fields/flat/polynomial_field.hpp"
#include "../../../geometry/planar/domain.hpp"
#include "../../../testing/flat_scatter.hpp"
#include "../../../testing/macros.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

using namespace geometry_art;
using namespace geometry_art::voronoi;
using namespace geometry_art::voronoi::flat;
using geometry::planar::Domain;
using geometry_art::testing::diagram_of;
using geometry_art::testing::flat_scatter;
using fields::flat::PolynomialField;

namespace {

constexpr double DISPLACEMENT = 1e-6;

std::vector<Vector3> cvt_gradient(const PolynomialField& field, const std::vector<Vector2>& sites) {
    CapacityConstrainedLagrangian<PolynomialField> lagrangian(
        field,
        field.total_mass() / static_cast<double>(sites.size())
    );

    return lagrangian
        .evaluate(*diagram_of(sites, field.domain()), std::vector<double>(sites.size(), 0.0), 0.0)
        .site_gradients;
}

void expect_matches_finite_differences(const PolynomialField& field, size_t count) {
    const Domain& domain = field.domain();
    std::vector<Vector2> sites = flat_scatter(count, domain.width, domain.height);
    auto diagram = diagram_of(sites, domain);

    HessianBlocks blocks = CvtHessian<PolynomialField>(field).assemble(*diagram);

    for (int trial = 0; trial < 3; ++trial) {
        std::vector<Vector3> direction(count);
        std::vector<Vector2> forward = sites;
        std::vector<Vector2> backward = sites;

        for (size_t k = 0; k < count; ++k) {
            direction[k] = Vector3(
                std::sin(1.3 * static_cast<double>(k) + trial),
                std::cos(2.1 * static_cast<double>(k) - trial),
                0.0
            );
            forward[k] += DISPLACEMENT * Vector2(direction[k].x(), direction[k].y());
            backward[k] -= DISPLACEMENT * Vector2(direction[k].x(), direction[k].y());
        }

        std::vector<Vector3> analytic = blocks.multiply(direction);
        std::vector<Vector3> ahead = cvt_gradient(field, forward);
        std::vector<Vector3> behind = cvt_gradient(field, backward);

        double scale = 0.0;

        for (size_t k = 0; k < count; ++k) {
            scale = std::max(scale, ((ahead[k] - behind[k]) / (2.0 * DISPLACEMENT)).norm());
        }

        for (size_t k = 0; k < count; ++k) {
            Vector3 difference = (ahead[k] - behind[k]) / (2.0 * DISPLACEMENT);

            EXPECT_LT((analytic[k] - difference).norm(), 1e-4 * scale) <<
                "site " << k << " trial " << trial;
        }
    }
}

} // namespace

TEST(FlatCvtHessianTest, MatchesFiniteDifferences) {
    expect_matches_finite_differences(PolynomialField::constant(1.0, Domain::torus(1.0, 1.0)), 12);
}

TEST(FlatCvtHessianTest, MatchesFiniteDifferencesOnAnElongatedTorus) {
    expect_matches_finite_differences(PolynomialField::constant(1.0, Domain::torus(3.0, 0.5)), 8);
}

TEST(FlatCvtHessianTest, IsSymmetric) {
    Domain domain = Domain::torus(1.0, 1.0);
    PolynomialField field = PolynomialField::constant(1.0, domain);
    std::vector<Vector2> sites = flat_scatter(13, 1.0, 1.0);
    auto diagram = diagram_of(sites, domain);

    HessianBlocks blocks = CvtHessian<PolynomialField>(field).assemble(*diagram);

    std::vector<Vector3> first(sites.size());
    std::vector<Vector3> second(sites.size());

    for (size_t k = 0; k < sites.size(); ++k) {
        first[k] = Vector3(std::sin(1.1 * static_cast<double>(k)), std::cos(0.6 * static_cast<double>(k)), 0.0);
        second[k] = Vector3(0.2, std::sin(2.3 * static_cast<double>(k)), 0.0);
    }

    std::vector<Vector3> applied_first = blocks.multiply(first);
    std::vector<Vector3> applied_second = blocks.multiply(second);

    double forward = 0.0;
    double backward = 0.0;

    for (size_t k = 0; k < sites.size(); ++k) {
        forward += second[k].dot(applied_first[k]);
        backward += first[k].dot(applied_second[k]);
    }

    EXPECT_NEAR(forward, backward, 1e-10);
}

// The energy's curvature reads only the bisectors; a wall neither sweeps
// nor contributes a relative second moment.
TEST(PlaneCvtHessianTest, MatchesFiniteDifferences) {
    expect_matches_finite_differences(PolynomialField::constant(1.0, Domain::plane(1.0, 1.0)), 12);
}

TEST(PlaneCvtHessianTest, MatchesFiniteDifferencesOnAGradient) {
    expect_matches_finite_differences(PolynomialField::linear(1.0, Vector2(1.0, 0.0), Domain::plane(2.0, 1.0)), 10);
}
