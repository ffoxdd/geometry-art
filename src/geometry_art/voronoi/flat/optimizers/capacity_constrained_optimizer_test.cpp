#include "capacity_constrained_optimizer.hpp"
#include "lloyd_optimizer.hpp"
#include "../core/torus.hpp"
#include "../../../fields/flat/constant_field.hpp"
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
using geometry_art::testing::flat_scatter;
using geometry_art::testing::scattered_torus;
using geometry_art::testing::torus_of;
using fields::flat::ConstantField;

namespace {

CapacityConstrainedParameters newton_parameters() {
    CapacityConstrainedParameters parameters;
    parameters.inner_solver = "newton";
    return parameters;
}

} // namespace

TEST(FlatOptimizerTest, LloydReducesTheCentroidDeviation) {
    ConstantField field(1.0, 1.0, 1.0);
    LloydOptimizer<ConstantField> lloyd(scattered_torus(12, 1.0, 1.0), field, 5, noop_callback());

    auto torus = lloyd.optimize();

    EXPECT_EQ(torus->size(), 12u);
    EXPECT_LT(lloyd.final_deviation(), 0.05);
}

TEST(FlatOptimizerTest, EXPENSIVE_EqualizesCapacitiesOnTheUnitTorus) {
    REQUIRE_EXPENSIVE();

    ConstantField field(1.0, 1.0, 1.0);
    CapacityConstrainedParameters parameters = newton_parameters();
    parameters.relative_capacity_tolerance = 1e-7;

    LloydOptimizer<ConstantField> lloyd(scattered_torus(20, 1.0, 1.0), field, 5, noop_callback());
    CapacityConstrainedOptimizer<ConstantField> optimizer(
        lloyd.optimize(),
        field,
        parameters,
        noop_callback()
    );

    auto torus = optimizer.optimize();
    const CapacityConstrainedReport& report = optimizer.report();

    EXPECT_TRUE(report.converged);
    EXPECT_LT(report.relative_rms_capacity_error, 1e-7);

    double target = field.total_mass() / 20.0;

    for (size_t index = 0; index < torus->size(); ++index) {
        EXPECT_NEAR(field.integrals(torus->cell(index)).mass, target, 1e-6 * target);
    }
}

TEST(FlatOptimizerTest, EXPENSIVE_EqualizesCapacitiesOnAnElongatedTorus) {
    REQUIRE_EXPENSIVE();

    ConstantField field(1.0, 3.0, 0.5);
    CapacityConstrainedParameters parameters = newton_parameters();
    parameters.relative_capacity_tolerance = 1e-7;

    LloydOptimizer<ConstantField> lloyd(scattered_torus(14, 3.0, 0.5), field, 5, noop_callback());
    CapacityConstrainedOptimizer<ConstantField> optimizer(
        lloyd.optimize(),
        field,
        parameters,
        noop_callback()
    );

    auto torus = optimizer.optimize();

    EXPECT_TRUE(optimizer.report().converged);
    EXPECT_LT(optimizer.report().relative_rms_capacity_error, 1e-7);
}
