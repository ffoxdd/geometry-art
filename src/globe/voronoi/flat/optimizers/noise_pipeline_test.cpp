#include "capacity_constrained_optimizer.hpp"
#include "lloyd_optimizer.hpp"
#include "../core/torus.hpp"
#include "../../../fields/flat/noise_field.hpp"
#include "../../../fields/flat/piecewise_polynomial_field.hpp"
#include "../../../math/interval.hpp"
#include "../../../testing/macros.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <cstddef>
#include <memory>

using namespace globe;
using namespace globe::voronoi;
using namespace globe::voronoi::flat;
using fields::flat::NoiseField;
using fields::flat::PiecewisePolynomialField;
using globe::math::Interval;

namespace {

std::unique_ptr<Torus> scattered_torus(size_t count, double width, double height) {
    auto torus = std::make_unique<Torus>(width, height);
    double golden = 0.6180339887498949;

    for (size_t k = 0; k < count; ++k) {
        double x = std::fmod(0.13 + golden * static_cast<double>(k), 1.0) * width;
        double y = (static_cast<double>(k) + 0.5) / static_cast<double>(count) * height;
        torus->insert(Vector2(x, y));
    }

    return torus;
}

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
