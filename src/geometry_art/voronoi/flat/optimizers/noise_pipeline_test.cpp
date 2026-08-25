#include "capacity_constrained_optimizer.hpp"
#include "lloyd_optimizer.hpp"
#include "../core/torus.hpp"
#include "../../../fields/flat/noise_field.hpp"
#include "../../../fields/flat/piecewise_polynomial_field.hpp"
#include "../../../math/interval.hpp"
#include "../../../testing/flat_scatter.hpp"
#include "../../../testing/macros.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <cstddef>
#include <memory>

using namespace geometry_art;
using namespace geometry_art::voronoi;
using namespace geometry_art::voronoi::flat;
using geometry_art::testing::flat_scatter;
using geometry_art::testing::scattered_torus;
using geometry_art::testing::torus_of;
using fields::flat::NoiseField;
using fields::flat::PiecewisePolynomialField;
using geometry_art::math::Interval;

namespace {

} // namespace

// The full flat pipeline on a real density: periodic noise sampled onto the
// piecewise representation, Lloyd warm start, then the capacity solve.
TEST(FlatNoisePipelineTest, EXPENSIVE_EqualizesCapacitiesUnderNoise) {
    REQUIRE_EXPENSIVE();

    NoiseField noise(2.0, 1.0, Interval(0.2, 1.0));
    auto field = PiecewisePolynomialField::sample(2.0, 1.0, 24, 12, noise);

    CapacityConstrainedParameters parameters;
    parameters.inner_solver = "newton";
    parameters.relative_capacity_tolerance = 1e-7;

    LloydOptimizer<PiecewisePolynomialField> lloyd(scattered_torus(24, 2.0, 1.0), field, 5, noop_callback());
    CapacityConstrainedOptimizer<PiecewisePolynomialField> optimizer(
        lloyd.optimize(),
        field,
        parameters,
        noop_callback()
    );

    auto torus = optimizer.optimize();

    EXPECT_TRUE(optimizer.report().converged);
    EXPECT_LT(optimizer.report().relative_rms_capacity_error, 1e-7);

    double target = field.total_mass() / 24.0;

    for (size_t index = 0; index < torus->size(); ++index) {
        EXPECT_NEAR(field.integrals(torus->cell(index)).mass, target, 1e-5 * target);
    }
}
