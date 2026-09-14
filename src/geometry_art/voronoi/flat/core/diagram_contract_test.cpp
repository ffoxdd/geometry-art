#include "diagram.hpp"
#include "../optimizers/capacity_constrained_lagrangian.hpp"
#include "../../../fields/flat/polynomial_field.hpp"
#include "../../../geometry/planar/domain.hpp"
#include "../../../testing/contracts/diagram_contract.hpp"
#include "../../../testing/contracts/state_contract.hpp"
#include "../../../testing/flat_scatter.hpp"
#include "../../../testing/macros.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

using namespace geometry_art;
using namespace geometry_art::testing::contracts;
using fields::flat::PolynomialField;
using geometry::planar::Domain;
using voronoi::flat::CapacityConstrainedLagrangian;
using voronoi::flat::Diagram;

namespace {

constexpr double WIDTH = 1.5;
constexpr double HEIGHT = 1.0;

template<typename DomainChoice>
struct FlatDiagramTraits {
    using DiagramType = Diagram;

    static std::unique_ptr<Diagram> scattered(size_t count) {
        return geometry_art::testing::scattered_diagram(count, DomainChoice::domain());
    }

    static double domain_area() { return WIDTH * HEIGHT; }
    static double cell_area(const Diagram& diagram, size_t index) { return diagram.cell(index).area(); }
    static Vector3 site(const Diagram& diagram, size_t index) { return diagram.site_vector(index); }

    static std::vector<EdgeRecord> edges(const Diagram& diagram, size_t index) {
        std::vector<EdgeRecord> records;

        for (const auto& edge : diagram.cell_edges(index)) {
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

template<typename DomainChoice>
struct FlatStateTraits {
    static voronoi::DiagramState scattered_state(size_t count) {
        PolynomialField field = PolynomialField::constant(0.8, DomainChoice::domain());
        CapacityConstrainedLagrangian<PolynomialField> lagrangian(field, field.total_mass() / count);
        return lagrangian.diagram_state(*geometry_art::testing::scattered_diagram(count, DomainChoice::domain()));
    }

    static double total_mass() { return 0.8 * WIDTH * HEIGHT; }
};

struct TorusChoice {
    static Domain domain() { return Domain::torus(WIDTH, HEIGHT); }
};

struct CylinderChoice {
    static Domain domain() { return Domain::cylinder(WIDTH, HEIGHT); }
};

struct PlaneChoice {
    static Domain domain() { return Domain::plane(WIDTH, HEIGHT); }
};

} // namespace

INSTANTIATE_TYPED_TEST_SUITE_P(Torus, DiagramContract, FlatDiagramTraits<TorusChoice>);
INSTANTIATE_TYPED_TEST_SUITE_P(Torus, StateContract, FlatStateTraits<TorusChoice>);
INSTANTIATE_TYPED_TEST_SUITE_P(Cylinder, DiagramContract, FlatDiagramTraits<CylinderChoice>);
INSTANTIATE_TYPED_TEST_SUITE_P(Cylinder, StateContract, FlatStateTraits<CylinderChoice>);
INSTANTIATE_TYPED_TEST_SUITE_P(Plane, DiagramContract, FlatDiagramTraits<PlaneChoice>);
INSTANTIATE_TYPED_TEST_SUITE_P(Plane, StateContract, FlatStateTraits<PlaneChoice>);
