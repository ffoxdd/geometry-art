#include "torus.hpp"
#include "../optimizers/capacity_constrained_lagrangian.hpp"
#include "../../../fields/flat/constant_field.hpp"
#include "../../../testing/contracts/diagram_contract.hpp"
#include "../../../testing/contracts/state_contract.hpp"
#include "../../../testing/macros.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

using namespace globe;
using namespace globe::testing::contracts;
using fields::flat::ConstantField;
using voronoi::flat::CapacityConstrainedLagrangian;
using voronoi::flat::Torus;

namespace {

constexpr double WIDTH = 1.5;
constexpr double HEIGHT = 1.0;

std::unique_ptr<Torus> scattered_torus(size_t count) {
    auto torus = std::make_unique<Torus>(WIDTH, HEIGHT);

    for (size_t k = 0; k < count; ++k) {
        double x = std::fmod(0.13 + 0.7548776662466927 * static_cast<double>(k), 1.0) * WIDTH;
        double y = std::fmod(0.41 + 0.5698402909980532 * static_cast<double>(k), 1.0) * HEIGHT;
        torus->insert(Vector2(x, y));
    }

    return torus;
}

struct TorusDiagramTraits {
    using DiagramType = Torus;

    static std::unique_ptr<Torus> scattered(size_t count) { return scattered_torus(count); }
    static double domain_area() { return WIDTH * HEIGHT; }
    static double cell_area(const Torus& torus, size_t index) { return torus.cell(index).area(); }
    static Vector3 site(const Torus& torus, size_t index) { return torus.site_vector(index); }

    static std::vector<EdgeRecord> edges(const Torus& torus, size_t index) {
        std::vector<EdgeRecord> records;

        for (const auto& edge : torus.cell_edges(index)) {
            records.push_back(EdgeRecord{
                edge.neighbor_index,
                Vector3(edge.neighbor_position.x(), edge.neighbor_position.y(), 0.0),
                Vector3(edge.boundary.source().x(), edge.boundary.source().y(), 0.0),
                Vector3(edge.boundary.target().x(), edge.boundary.target().y(), 0.0)
            });
        }

        return records;
    }
};

struct TorusStateTraits {
    static voronoi::DiagramState scattered_state(size_t count) {
        ConstantField field(0.8, WIDTH, HEIGHT);
        CapacityConstrainedLagrangian<ConstantField> lagrangian(field, field.total_mass() / count);
        return lagrangian.diagram_state(*scattered_torus(count));
    }

    static double total_mass() { return 0.8 * WIDTH * HEIGHT; }
};

} // namespace

INSTANTIATE_TYPED_TEST_SUITE_P(Torus, DiagramContract, TorusDiagramTraits);
INSTANTIATE_TYPED_TEST_SUITE_P(Torus, StateContract, TorusStateTraits);
