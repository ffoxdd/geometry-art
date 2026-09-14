#include "capacity_constrained_optimizer.hpp"
#include "lloyd_optimizer.hpp"
#include "../core/diagram.hpp"
#include "../../../fields/flat/noise_field.hpp"
#include "../../../fields/flat/piecewise_polynomial_field.hpp"
#include "../../../geometry/planar/domain.hpp"
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
using geometry::planar::Domain;
using geometry_art::testing::scattered_diagram;
using fields::flat::NoiseField;
using fields::flat::PiecewisePolynomialField;
using geometry_art::math::Interval;

namespace {

// The full flat pipeline on a real density: noise sampled onto the
// piecewise representation, Lloyd warm start, then the capacity solve.
void expect_equalizes_capacities_under_noise(const Domain& domain) {
    NoiseField noise(domain.width, domain.height, Interval(0.2, 1.0));
    auto field = PiecewisePolynomialField::sample(domain, 24, 12, noise);

    CapacityConstrainedParameters parameters;
    parameters.inner_solver = "newton";
    parameters.relative_capacity_tolerance = 1e-7;

    LloydOptimizer<PiecewisePolynomialField> lloyd(scattered_diagram(24, domain), field, 5, noop_callback());
    CapacityConstrainedOptimizer<PiecewisePolynomialField> optimizer(
        lloyd.optimize(),
        field,
        parameters,
        noop_callback()
    );

    auto diagram = optimizer.optimize();

    EXPECT_TRUE(optimizer.report().converged);
    EXPECT_LT(optimizer.report().relative_rms_capacity_error, 1e-7);

    double target = field.total_mass() / 24.0;

    for (size_t index = 0; index < diagram->size(); ++index) {
        EXPECT_NEAR(field.integrals(diagram->cell(index)).mass, target, 1e-5 * target);
    }
}

} // namespace

TEST(FlatNoisePipelineTest, EXPENSIVE_EqualizesCapacitiesUnderNoise) {
    REQUIRE_EXPENSIVE();
    expect_equalizes_capacities_under_noise(Domain::torus(2.0, 1.0));
}

TEST(PlaneNoisePipelineTest, EXPENSIVE_EqualizesCapacitiesUnderNoise) {
    REQUIRE_EXPENSIVE();
    expect_equalizes_capacities_under_noise(Domain::plane(2.0, 1.0));
}
