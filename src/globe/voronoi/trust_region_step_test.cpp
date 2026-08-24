#include "trust_region_step.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <vector>

using namespace globe;
using namespace globe::voronoi;

namespace {

// A diagonal Hessian makes the exact minimiser available in closed form,
// so the solver can be checked without any geometry in the way.
class DiagonalOperator {
 public:
    explicit DiagonalOperator(std::vector<double> weights) : _weights(std::move(weights)) {}

    [[nodiscard]] std::vector<Vector3> multiply(const std::vector<Vector3>& directions) const {
        std::vector<Vector3> result(directions.size());

        for (size_t i = 0; i < directions.size(); ++i) {
            result[i] = _weights[i] * directions[i];
        }

        return result;
    }

 private:
    std::vector<double> _weights;
};

double norm(const std::vector<Vector3>& value) {
    double sum = 0.0;

    for (const Vector3& entry : value) {
        sum += entry.squaredNorm();
    }

    return std::sqrt(sum);
}

} // namespace

TEST(TrustRegionStepTest, ReachesTheExactMinimiserInsideTheRegion) {
    std::vector<Vector3> gradient{Vector3(1.0, -2.0, 0.5), Vector3(0.25, 0.0, -1.5)};
    DiagonalOperator hessian({2.0, 4.0});

    TrustRegionStep::Result result = TrustRegionStep(50, 1e-12).solve(gradient, hessian, 100.0);

    EXPECT_FALSE(result.hit_boundary);
    EXPECT_NEAR((result.step[0] + gradient[0] / 2.0).norm(), 0.0, 1e-10);
    EXPECT_NEAR((result.step[1] + gradient[1] / 4.0).norm(), 0.0, 1e-10);
    EXPECT_GT(result.predicted_decrease, 0.0);
}

TEST(TrustRegionStepTest, StopsOnTheBoundaryWhenTheRegionIsSmall) {
    std::vector<Vector3> gradient{Vector3(1.0, -2.0, 0.5), Vector3(0.25, 0.0, -1.5)};
    DiagonalOperator hessian({2.0, 4.0});
    double radius = 0.1;

    TrustRegionStep::Result result = TrustRegionStep(50, 1e-12).solve(gradient, hessian, radius);

    EXPECT_TRUE(result.hit_boundary);
    EXPECT_NEAR(norm(result.step), radius, 1e-12);
    EXPECT_GT(result.predicted_decrease, 0.0);
}

TEST(TrustRegionStepTest, FollowsNegativeCurvatureToTheBoundary) {
    std::vector<Vector3> gradient{Vector3(0.3, 0.0, 0.0)};
    DiagonalOperator hessian({-1.0});
    double radius = 2.0;

    TrustRegionStep::Result result = TrustRegionStep(50, 1e-12).solve(gradient, hessian, radius);

    EXPECT_TRUE(result.hit_boundary);
    EXPECT_NEAR(norm(result.step), radius, 1e-12);
    EXPECT_LT(result.step[0].x(), 0.0);
}

TEST(TrustRegionStepTest, TakesNoStepWhenTheGradientVanishes) {
    std::vector<Vector3> gradient{Vector3::Zero(), Vector3::Zero()};
    DiagonalOperator hessian({2.0, 4.0});

    TrustRegionStep::Result result = TrustRegionStep(50, 1e-12).solve(gradient, hessian, 1.0);

    EXPECT_EQ(result.iterations, 0u);
    EXPECT_NEAR(norm(result.step), 0.0, 1e-15);
}
