#include "gtest/gtest.h"
#include "harmonic_hit.hpp"
#include "../../math/interval.hpp"
#include "../../math/interval_sampler/uniform_interval_sampler.hpp"
#include "../../testing/macros.hpp"
#include <algorithm>
#include <cmath>
#include <vector>

using namespace geometry_art::dla;
using geometry_art::Interval;
using geometry_art::UniformIntervalSampler;

TEST(HarmonicHitTest, ReturnProbabilityIsRadiusOverDistance) {
    EXPECT_DOUBLE_EQ(return_probability(4.0, 2.0), 0.5);
    EXPECT_DOUBLE_EQ(return_probability(10.0, 1.0), 0.1);
}

TEST(HarmonicHitTest, ZeroDrawLandsOnTheFarSide) {
    EXPECT_NEAR(harmonic_hit_cosine(2.0, 1.0, 0.0), -1.0, 1e-12);
}

TEST(HarmonicHitTest, FullDrawLandsOnTheNearSide) {
    EXPECT_NEAR(harmonic_hit_cosine(2.0, 1.0, 1.0), 1.0, 1e-12);
}

TEST(HarmonicHitTest, MedianDrawMatchesTheClosedForm) {
    EXPECT_NEAR(harmonic_hit_cosine(2.0, 1.0, 0.5), 0.6875, 1e-12);
}

TEST(HarmonicHitTest, CosineGrowsWithTheDraw) {
    EXPECT_LT(harmonic_hit_cosine(3.0, 1.0, 0.25), harmonic_hit_cosine(3.0, 1.0, 0.75));
}

namespace {

// The polar CDF of the exterior Poisson kernel, written out independently of
// the inverse the sampler uses.
double analytic_cdf(double distance, double radius, double cosine) {
    double axis_product = 2.0 * distance * radius;
    double squared_sum = distance * distance + radius * radius;

    return (distance * distance - radius * radius) / (2.0 * radius) * (
        1.0 / std::sqrt(squared_sum - axis_product * cosine) -
        1.0 / (distance + radius)
    );
}

} // namespace

TEST(HarmonicHitTest, EXPENSIVE_SampledCosinesMatchTheAnalyticDistribution) {
    REQUIRE_EXPENSIVE();

    constexpr double distance = 2.0;
    constexpr double radius = 1.0;
    constexpr size_t sample_count = 100000;

    UniformIntervalSampler sampler(42);
    std::vector<double> cosines;
    cosines.reserve(sample_count);

    for (size_t i = 0; i < sample_count; ++i) {
        cosines.push_back(harmonic_hit_cosine(distance, radius, sampler.sample(geometry_art::UNIT_INTERVAL)));
    }

    for (double threshold : {-0.5, 0.0, 0.5, 0.9}) {
        double below = static_cast<double>(std::count_if(
            cosines.begin(), cosines.end(),
            [threshold](double cosine) { return cosine <= threshold; }
        ));

        EXPECT_NEAR(below / sample_count, analytic_cdf(distance, radius, threshold), 0.01)
            << "empirical CDF diverges from the Poisson kernel at cosine " << threshold;
    }
}
