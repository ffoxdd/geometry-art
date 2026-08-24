#include "capacity_constrained_lagrangian.hpp"
#include "../core/torus.hpp"
#include "../../../fields/flat/constant_field.hpp"
#include "../../../testing/macros.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

using namespace globe;
using namespace globe::voronoi;
using namespace globe::voronoi::flat;
using fields::flat::ConstantField;

namespace {

constexpr double DISPLACEMENT = 1e-6;

std::vector<Vector2> scattered_sites(size_t count, double width, double height) {
    std::vector<Vector2> sites;
    double golden = 0.6180339887498949;

    for (size_t k = 0; k < count; ++k) {
        double x = std::fmod(0.13 + golden * static_cast<double>(k), 1.0) * width;
        double y = (static_cast<double>(k) + 0.5) / static_cast<double>(count) * height;
        sites.emplace_back(x, y);
    }

    return sites;
}

std::unique_ptr<Torus> torus_of(const std::vector<Vector2>& sites, double width, double height) {
    auto torus = std::make_unique<Torus>(width, height);

    for (const Vector2& site : sites) {
        torus->insert(site);
    }

    return torus;
}

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
void expect_gradient_matches_differences(size_t count, double width, double height) {
    ConstantField field(1.0, width, height);
    std::vector<Vector2> sites = scattered_sites(count, width, height);
    auto torus = torus_of(sites, width, height);

    CapacityConstrainedLagrangian<ConstantField> lagrangian(
        field,
        field.total_mass() / static_cast<double>(count)
    );

    std::vector<double> multipliers = arbitrary_multipliers(count);
    double penalty = 0.7;
    LagrangianEvaluation evaluation = lagrangian.evaluate(*torus, multipliers, penalty);

    for (size_t k = 0; k < count; ++k) {
        for (int axis = 0; axis < 2; ++axis) {
            Vector2 step = Vector2::Zero();
            step[axis] = DISPLACEMENT;

            std::vector<Vector2> forward = sites;
            std::vector<Vector2> backward = sites;
            forward[k] += step;
            backward[k] -= step;

            double ahead = lagrangian
                .evaluate(*torus_of(forward, width, height), multipliers, penalty)
                .value;
            double behind = lagrangian
                .evaluate(*torus_of(backward, width, height), multipliers, penalty)
                .value;

            double difference = (ahead - behind) / (2.0 * DISPLACEMENT);

            EXPECT_NEAR(evaluation.site_gradients[k][axis], difference, 1e-5) <<
                "site " << k << " axis " << axis;
        }
    }
}

} // namespace

TEST(FlatLagrangianTest, CapacitiesPartitionTheTotalMass) {
    ConstantField field(2.0, 1.5, 1.0);
    auto torus = torus_of(scattered_sites(12, 1.5, 1.0), 1.5, 1.0);
    CapacityConstrainedLagrangian<ConstantField> lagrangian(field, field.total_mass() / 12.0);

    DiagramState state = lagrangian.diagram_state(*torus);
    double total = 0.0;

    for (const CellState& cell : state.cells) {
        total += cell.mass;
    }

    EXPECT_NEAR(total, field.total_mass(), 1e-9);
}

TEST(FlatLagrangianTest, RelativeMomentsSatisfyTheSeparationIdentity) {
    ConstantField field(1.0, 1.0, 1.0);
    auto torus = torus_of(scattered_sites(9, 1.0, 1.0), 1.0, 1.0);
    CapacityConstrainedLagrangian<ConstantField> lagrangian(field, field.total_mass() / 9.0);

    DiagramState state = lagrangian.diagram_state(*torus);

    for (size_t k = 0; k < state.edges.size(); ++k) {
        for (const EdgeState& edge : state.edges[k]) {
            Vector3 difference = edge.moment_about_own - edge.moment_about_neighbor;

            EXPECT_NEAR(difference.norm(), edge.separation * edge.mass, 1e-9);
        }
    }
}

TEST(FlatLagrangianTest, EnergyOfAGridIsTheClosedFormValue) {
    // A 2x2 grid on the unit torus: each cell is a 0.5 square centered on
    // its site, whose integral of |x - s|^2 is (side^4)/6.
    ConstantField field(1.0, 1.0, 1.0);
    Torus torus(1.0, 1.0);

    for (int row = 0; row < 2; ++row) {
        for (int column = 0; column < 2; ++column) {
            torus.insert(Vector2(0.25 + 0.5 * column, 0.25 + 0.5 * row));
        }
    }

    CapacityConstrainedLagrangian<ConstantField> lagrangian(field, 0.25);
    LagrangianEvaluation evaluation = lagrangian.evaluate(torus, std::vector<double>(4, 0.0), 0.0);

    EXPECT_NEAR(evaluation.cvt_energy, 4.0 * std::pow(0.5, 4) / 6.0, 1e-9);
    EXPECT_NEAR(evaluation.root_mean_square_capacity_error(), 0.0, 1e-9);
}

TEST(FlatLagrangianTest, GradientMatchesCentralDifferences) {
    expect_gradient_matches_differences(10, 1.0, 1.0);
}

// An elongated torus forces bisector pairs that differ only by the period
// offset, so the slot keys' period component is what this exercises.
TEST(FlatLagrangianTest, GradientMatchesCentralDifferencesOnAnElongatedTorus) {
    expect_gradient_matches_differences(6, 3.0, 0.5);
}
