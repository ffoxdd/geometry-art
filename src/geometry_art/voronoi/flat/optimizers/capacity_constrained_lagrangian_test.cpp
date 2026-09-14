#include "capacity_constrained_lagrangian.hpp"
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
using geometry_art::testing::diagram_of;
using geometry_art::testing::flat_scatter;
using fields::flat::PolynomialField;

namespace {

constexpr double DISPLACEMENT = 1e-6;

std::vector<double> arbitrary_multipliers(size_t count) {
    std::vector<double> multipliers(count);

    for (size_t k = 0; k < count; ++k) {
        multipliers[k] = std::sin(2.3 * static_cast<double>(k) + 0.4);
    }

    return multipliers;
}

// The value is a chart-invariant scalar of the sites, so its gradient is
// checked against central differences of the value itself, with the diagram
// rebuilt at every displaced configuration.
void expect_gradient_matches_differences(const PolynomialField& field, size_t count) {
    const Domain& domain = field.domain();
    std::vector<Vector2> sites = flat_scatter(count, domain.width, domain.height);
    auto diagram = diagram_of(sites, domain);

    CapacityConstrainedLagrangian<PolynomialField> lagrangian(
        field,
        field.total_mass() / static_cast<double>(count)
    );

    std::vector<double> multipliers = arbitrary_multipliers(count);
    double penalty = 0.7;
    LagrangianEvaluation evaluation = lagrangian.evaluate(*diagram, multipliers, penalty);

    for (size_t k = 0; k < count; ++k) {
        for (int axis = 0; axis < 2; ++axis) {
            Vector2 step = Vector2::Zero();
            step[axis] = DISPLACEMENT;

            std::vector<Vector2> forward = sites;
            std::vector<Vector2> backward = sites;
            forward[k] += step;
            backward[k] -= step;

            double ahead = lagrangian
                .evaluate(*diagram_of(forward, domain), multipliers, penalty)
                .value;
            double behind = lagrangian
                .evaluate(*diagram_of(backward, domain), multipliers, penalty)
                .value;

            double difference = (ahead - behind) / (2.0 * DISPLACEMENT);

            EXPECT_NEAR(evaluation.site_gradients[k][axis], difference, 1e-5) <<
                "site " << k << " axis " << axis;
        }
    }
}

std::vector<Vector2> grid_sites() {
    std::vector<Vector2> sites;

    for (int row = 0; row < 2; ++row) {
        for (int column = 0; column < 2; ++column) {
            sites.emplace_back(0.25 + 0.5 * column, 0.25 + 0.5 * row);
        }
    }

    return sites;
}

} // namespace

TEST(FlatLagrangianTest, CapacitiesPartitionTheTotalMass) {
    Domain domain = Domain::torus(1.5, 1.0);
    PolynomialField field = PolynomialField::constant(2.0, domain);
    auto diagram = diagram_of(flat_scatter(12, 1.5, 1.0), domain);
    CapacityConstrainedLagrangian<PolynomialField> lagrangian(field, field.total_mass() / 12.0);

    DiagramState state = lagrangian.diagram_state(*diagram);
    double total = 0.0;

    for (const CellState& cell : state.cells) {
        total += cell.mass;
    }

    EXPECT_NEAR(total, field.total_mass(), 1e-9);
}

TEST(FlatLagrangianTest, RelativeMomentsSatisfyTheSeparationIdentity) {
    Domain domain = Domain::torus(1.0, 1.0);
    PolynomialField field = PolynomialField::constant(1.0, domain);
    auto diagram = diagram_of(flat_scatter(9, 1.0, 1.0), domain);
    CapacityConstrainedLagrangian<PolynomialField> lagrangian(field, field.total_mass() / 9.0);

    DiagramState state = lagrangian.diagram_state(*diagram);

    for (size_t k = 0; k < state.edges.size(); ++k) {
        for (const EdgeState& edge : state.edges[k]) {
            Vector3 difference = edge.moment_about_own - edge.moment_about_neighbor;

            EXPECT_NEAR(difference.norm(), edge.separation * edge.mass, 1e-9);
        }
    }
}

// A 2x2 grid on the unit torus: each cell is a 0.5 square centered on its
// site, whose integral of |x - s|^2 is (side^4)/6.
TEST(FlatLagrangianTest, EnergyOfAGridIsTheClosedFormValue) {
    Domain domain = Domain::torus(1.0, 1.0);
    PolynomialField field = PolynomialField::constant(1.0, domain);
    Diagram diagram(domain, grid_sites());

    CapacityConstrainedLagrangian<PolynomialField> lagrangian(field, 0.25);
    LagrangianEvaluation evaluation = lagrangian.evaluate(diagram, std::vector<double>(4, 0.0), 0.0);

    EXPECT_NEAR(evaluation.cvt_energy, 4.0 * std::pow(0.5, 4) / 6.0, 1e-9);
    EXPECT_NEAR(evaluation.root_mean_square_capacity_error(), 0.0, 1e-9);
}

// The same grid on the plane has the same cells, and so the same energy.
TEST(PlaneLagrangianTest, EnergyOfAGridIsTheClosedFormValue) {
    Domain domain = Domain::plane(1.0, 1.0);
    PolynomialField field = PolynomialField::constant(1.0, domain);
    Diagram diagram(domain, grid_sites());

    CapacityConstrainedLagrangian<PolynomialField> lagrangian(field, 0.25);
    LagrangianEvaluation evaluation = lagrangian.evaluate(diagram, std::vector<double>(4, 0.0), 0.0);

    EXPECT_NEAR(evaluation.cvt_energy, 4.0 * std::pow(0.5, 4) / 6.0, 1e-9);
    EXPECT_NEAR(evaluation.root_mean_square_capacity_error(), 0.0, 1e-9);
}

TEST(FlatLagrangianTest, GradientMatchesCentralDifferences) {
    expect_gradient_matches_differences(PolynomialField::constant(1.0, Domain::torus(1.0, 1.0)), 10);
}

// An elongated torus forces bisector pairs that differ only by the period
// offset, so the slot keys' period component is what this exercises.
TEST(FlatLagrangianTest, GradientMatchesCentralDifferencesOnAnElongatedTorus) {
    expect_gradient_matches_differences(PolynomialField::constant(1.0, Domain::torus(3.0, 0.5)), 6);
}

// Wall edges do not sweep, so the gradient of a walled cell reads only its
// bisectors; the differences would expose any wall term wrongly counted.
TEST(PlaneLagrangianTest, GradientMatchesCentralDifferences) {
    expect_gradient_matches_differences(PolynomialField::constant(1.0, Domain::plane(1.0, 1.0)), 10);
}

TEST(PlaneLagrangianTest, GradientMatchesCentralDifferencesOnAGradient) {
    expect_gradient_matches_differences(PolynomialField::linear(1.0, Vector2(1.0, 0.0), Domain::plane(2.0, 1.0)), 8);
}
