#include "sphere.hpp"
#include "../optimizers/capacity_constrained_lagrangian.hpp"
#include "../../../fields/spherical/polynomial_field.hpp"
#include "../../../testing/contracts/diagram_contract.hpp"
#include "../../../testing/contracts/state_contract.hpp"
#include "../../../testing/macros.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

using namespace geometry_art;
using namespace geometry_art::testing::contracts;
using fields::spherical::PolynomialField;
using voronoi::spherical::CapacityConstrainedLagrangian;
using voronoi::spherical::Sphere;

namespace {

std::unique_ptr<Sphere> fibonacci_sphere(size_t count) {
    auto sphere = std::make_unique<Sphere>();
    double golden_angle = M_PI * (std::sqrt(5.0) - 1.0);

    for (size_t k = 0; k < count; ++k) {
        double y = 1.0 - 2.0 * (static_cast<double>(k) + 0.5) / static_cast<double>(count);
        double radius = std::sqrt(1.0 - y * y);
        double angle = golden_angle * static_cast<double>(k);
        sphere->insert(cgal::to_point(VectorS2(std::cos(angle) * radius, y, std::sin(angle) * radius)));
    }

    return sphere;
}

PolynomialField contract_field() {
    Eigen::Matrix3d quadratic = Eigen::Matrix3d::Zero();
    quadratic(2, 2) = -0.5;
    return PolynomialField::quadratic(1.0, Vector3(0.1, -0.2, 0.3), quadratic);
}

struct SphereDiagramTraits {
    using DiagramType = Sphere;

    static std::unique_ptr<Sphere> scattered(size_t count) { return fibonacci_sphere(count); }
    static double domain_area() { return 4.0 * M_PI; }
    static double cell_area(const Sphere& sphere, size_t index) { return sphere.cell(index).area(); }
    static Vector3 site(const Sphere& sphere, size_t index) { return to_vector3(sphere.site(index)); }

    static std::vector<EdgeRecord> edges(const Sphere& sphere, size_t index) {
        std::vector<EdgeRecord> records;

        for (const auto& edge : sphere.cell_edges(index)) {
            records.push_back(EdgeRecord{
                edge.neighbor_index,
                to_vector3(sphere.site(edge.neighbor_index)),
                Vector3(edge.arc.source()),
                Vector3(edge.arc.target())
            });
        }

        return records;
    }
};

struct SphereStateTraits {
    static voronoi::DiagramState scattered_state(size_t count) {
        PolynomialField field = contract_field();
        CapacityConstrainedLagrangian<PolynomialField> lagrangian(field, field.total_mass() / count);
        return lagrangian.sphere_state(*fibonacci_sphere(count));
    }

    static double total_mass() { return contract_field().total_mass(); }
};

} // namespace

INSTANTIATE_TYPED_TEST_SUITE_P(Sphere, DiagramContract, SphereDiagramTraits);
INSTANTIATE_TYPED_TEST_SUITE_P(Sphere, StateContract, SphereStateTraits);
