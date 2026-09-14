#include "capacity_constrained_optimizer.hpp"
#include "lloyd_optimizer.hpp"
#include "../core/diagram.hpp"
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
using geometry_art::testing::scattered_diagram;
using fields::flat::PolynomialField;

namespace {

CapacityConstrainedParameters newton_parameters() {
    CapacityConstrainedParameters parameters;
    parameters.inner_solver = "newton";
    parameters.relative_capacity_tolerance = 1e-7;
    return parameters;
}

void expect_equalizes_capacities(const PolynomialField& field, size_t count, double tolerance) {
    LloydOptimizer<PolynomialField> lloyd(scattered_diagram(count, field.domain()), field, 5, noop_callback());
    CapacityConstrainedOptimizer<PolynomialField> optimizer(
        lloyd.optimize(),
        field,
        newton_parameters(),
        noop_callback()
    );

    auto diagram = optimizer.optimize();
    const CapacityConstrainedReport& report = optimizer.report();

    EXPECT_TRUE(report.converged);
    EXPECT_LT(report.relative_rms_capacity_error, 1e-7);

    double target = field.total_mass() / static_cast<double>(count);

    for (size_t index = 0; index < diagram->size(); ++index) {
        EXPECT_NEAR(field.integrals(diagram->cell(index)).mass, target, tolerance * target);
    }
}

} // namespace

TEST(FlatOptimizerTest, LloydReducesTheCentroidDeviation) {
    Domain domain = Domain::torus(1.0, 1.0);
    PolynomialField field = PolynomialField::constant(1.0, domain);
    LloydOptimizer<PolynomialField> lloyd(scattered_diagram(12, domain), field, 5, noop_callback());

    auto diagram = lloyd.optimize();

    EXPECT_EQ(diagram->size(), 12u);
    EXPECT_LT(lloyd.final_deviation(), 0.05);
}

TEST(FlatOptimizerTest, EXPENSIVE_EqualizesCapacitiesOnTheUnitTorus) {
    REQUIRE_EXPENSIVE();
    expect_equalizes_capacities(PolynomialField::constant(1.0, Domain::torus(1.0, 1.0)), 20, 1e-6);
}

TEST(FlatOptimizerTest, EXPENSIVE_EqualizesCapacitiesOnAnElongatedTorus) {
    REQUIRE_EXPENSIVE();
    expect_equalizes_capacities(PolynomialField::constant(1.0, Domain::torus(3.0, 0.5)), 14, 1e-6);
}

TEST(PlaneOptimizerTest, LloydKeepsSitesInsideTheWalls) {
    Domain domain = Domain::plane(2.0, 1.0);
    PolynomialField field = PolynomialField::constant(1.0, domain);
    LloydOptimizer<PolynomialField> lloyd(scattered_diagram(12, domain), field, 5, noop_callback());

    auto diagram = lloyd.optimize();

    for (size_t index = 0; index < diagram->size(); ++index) {
        Vector2 site = diagram->site(index);
        EXPECT_GT(site.x(), 0.0);
        EXPECT_LT(site.x(), 2.0);
        EXPECT_GT(site.y(), 0.0);
        EXPECT_LT(site.y(), 1.0);
    }

    EXPECT_LT(lloyd.final_deviation(), 0.05);
}

TEST(PlaneOptimizerTest, EXPENSIVE_EqualizesCapacitiesOnAConstantDensity) {
    REQUIRE_EXPENSIVE();
    expect_equalizes_capacities(PolynomialField::constant(1.0, Domain::plane(2.0, 1.0)), 20, 1e-6);
}

TEST(CylinderOptimizerTest, EXPENSIVE_EqualizesCapacitiesOnAGradientUpTheHeight) {
    REQUIRE_EXPENSIVE();
    expect_equalizes_capacities(PolynomialField::linear(1.0, Vector2(0.0, 2.0), Domain::cylinder(2.0, 1.0)), 24, 1e-6);
}

// The gradient is the density the plane exists for: a torus cannot carry
// one, and every cell along the rise must hold the same mass.
TEST(PlaneOptimizerTest, EXPENSIVE_EqualizesCapacitiesOnAGradient) {
    REQUIRE_EXPENSIVE();
    expect_equalizes_capacities(PolynomialField::linear(1.0, Vector2(1.5, 0.0), Domain::plane(2.0, 1.0)), 24, 1e-6);
}
