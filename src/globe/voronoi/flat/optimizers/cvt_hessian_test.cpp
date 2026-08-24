#include "cvt_hessian.hpp"
#include "capacity_constrained_lagrangian.hpp"
#include "../core/torus.hpp"
#include "../../hessian_blocks.hpp"
#include "../../../fields/flat/constant_field.hpp"
#include "../../../testing/flat_scatter.hpp"
#include "../../../testing/macros.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

using namespace globe;
using namespace globe::voronoi;
using namespace globe::voronoi::flat;
using globe::testing::flat_scatter;
using globe::testing::scattered_torus;
using globe::testing::torus_of;
using fields::flat::ConstantField;

namespace {

constexpr double DISPLACEMENT = 1e-6;

std::vector<Vector3> cvt_gradient(
    const ConstantField& field,
    const std::vector<Vector2>& sites,
    double width,
    double height
) {
    CapacityConstrainedLagrangian<ConstantField> lagrangian(
        field,
        field.total_mass() / static_cast<double>(sites.size())
    );

    return lagrangian
        .evaluate(*torus_of(sites, width, height), std::vector<double>(sites.size(), 0.0), 0.0)
        .site_gradients;
}

void expect_matches_finite_differences(size_t count, double width, double height) {
    ConstantField field(1.0, width, height);
    std::vector<Vector2> sites = flat_scatter(count, width, height);
    auto torus = torus_of(sites, width, height);

    HessianBlocks blocks = CvtHessian<ConstantField>(field).assemble(*torus);

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
        std::vector<Vector3> ahead = cvt_gradient(field, forward, width, height);
        std::vector<Vector3> behind = cvt_gradient(field, backward, width, height);

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
    expect_matches_finite_differences(10, 1.0, 1.0);
}

TEST(FlatCvtHessianTest, MatchesFiniteDifferencesOnAnElongatedTorus) {
    expect_matches_finite_differences(6, 3.0, 0.5);
}

TEST(FlatCvtHessianTest, IsSymmetric) {
    ConstantField field(1.0, 1.0, 1.0);
    std::vector<Vector2> sites = flat_scatter(12, 1.0, 1.0);
    auto torus = torus_of(sites, 1.0, 1.0);

    HessianBlocks blocks = CvtHessian<ConstantField>(field).assemble(*torus);

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
    double scale = 0.0;

    for (size_t k = 0; k < sites.size(); ++k) {
        forward += second[k].dot(applied_first[k]);
        backward += first[k].dot(applied_second[k]);
        scale = std::max(scale, applied_first[k].norm());
    }

    EXPECT_NEAR(forward, backward, 1e-10 * std::max(scale, 1.0));
}
