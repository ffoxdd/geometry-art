#ifndef GEOMETRY_ART_VORONOI_FLAT_OPTIMIZERS_CAPACITY_HESSIAN_HPP_
#define GEOMETRY_ART_VORONOI_FLAT_OPTIMIZERS_CAPACITY_HESSIAN_HPP_

#include "../core/cut.hpp"
#include "../core/periodic_slots.hpp"
#include "../core/diagram.hpp"
#include "../../hessian_blocks.hpp"
#include "../../../fields/flat/field.hpp"
#include "../../../std_ext/parallel_for.hpp"
#include "../../../types.hpp"
#include <CGAL/assertions.h>
#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <cstddef>
#include <utility>
#include <variant>
#include <vector>

namespace geometry_art::voronoi::flat {

// The constraints' own curvature: second derivatives of the weighted cell
// masses. A bisector's sweep rate responds to site motion through three
// channels -- the line translates and rotates, which the density and its
// gradient feel along the edge; the separation in the denominator changes;
// and the two endpoints slide, each pinned by planar equidistance to one
// more site or by the wall it lies on. A self-edge's weight jump is zero,
// so self-edges contribute nothing and are skipped; a wall edge never
// sweeps and is not a bisector at all. An endpoint pinned by a wall has no
// fourth site, and its slot names the cell itself with an empty row.
template<fields::flat::Field FieldType>
class CapacityHessian {
 public:
    explicit CapacityHessian(FieldType field);

    [[nodiscard]] HessianBlocks assemble(const Diagram& diagram, const std::vector<double>& weights) const;

 private:
    struct EdgeContribution {
        std::array<size_t, 4> sites;
        std::array<Matrix3, 4> own_row;
        std::array<Matrix3, 4> neighbor_row;
        bool active;
    };

    // The linear response of the edge's density integrals -- its mass and
    // its raw first moment -- to each of the four sites' motions.
    struct Response {
        std::array<Eigen::RowVector3d, 4> mass;
        std::array<Matrix3, 4> moment;
    };

    FieldType _field;

    [[nodiscard]] EdgeContribution edge_contribution(const Diagram& diagram, size_t cell, const CellEdgeInfo& edge) const;

    void add_endpoint(
        Response& response,
        const Vector2& endpoint,
        const Vector2& own,
        const Vector2& neighbor,
        const Cut& cut,
        int opposite_slot,
        const Vector3& tangent,
        double sign
    ) const;

    [[nodiscard]] static size_t pinned_site(const Cut& cut, size_t cell);
    [[nodiscard]] static Matrix3 planar_identity();
    [[nodiscard]] static Vector3 planar(const Vector2& point);
};

template<fields::flat::Field FieldType>
CapacityHessian<FieldType>::CapacityHessian(FieldType field) :
    _field(std::move(field)) {
}

template<fields::flat::Field FieldType>
HessianBlocks CapacityHessian<FieldType>::assemble(
    const Diagram& diagram,
    const std::vector<double>& weights
) const {
    size_t count = diagram.size();
    CGAL_precondition(weights.size() == count);

    std::vector<std::vector<CellEdgeInfo>> cell_edges(count);

    for (size_t k = 0; k < count; ++k) {
        cell_edges[k] = diagram.cell_edges(k);
    }

    PeriodicSlots slots = PeriodicSlots::build(diagram, cell_edges);
    std::vector<EdgeContribution> contributions(slots.count());

    std_ext::parallel_for(slots.count(), [&](size_t slot) {
        auto [cell, position] = slots.representatives[slot];
        contributions[slot] = edge_contribution(diagram, cell, cell_edges[cell][position]);
    });

    HessianBlocks blocks;
    blocks.diagonal.assign(count, Matrix3::Zero());
    blocks.neighbors.resize(count);

    for (size_t k = 0; k < count; ++k) {
        blocks.neighbors[k].reserve(cell_edges[k].size());

        for (const CellEdgeInfo& edge : cell_edges[k]) {
            blocks.neighbors[k].push_back(NeighborBlock{edge.neighbor_index, Matrix3::Zero()});
        }
    }

    // Duplicate neighbor entries are fine as scatter targets: the operator
    // sums blocks linearly, so adding to the first entry is equivalent.
    auto scatter = [&](size_t row, size_t column, const Matrix3& block) {
        if (column == row) {
            blocks.diagonal[row] += block;
            return;
        }

        for (NeighborBlock& neighbor : blocks.neighbors[row]) {
            if (neighbor.neighbor_index == column) {
                neighbor.value += block;
                return;
            }
        }

        CGAL_assertion(false);
    };

    for (size_t slot = 0; slot < slots.count(); ++slot) {
        const EdgeContribution& contribution = contributions[slot];

        if (!contribution.active) {
            continue;
        }

        size_t own = contribution.sites[0];
        size_t neighbor = contribution.sites[1];
        double jump = weights[own] - weights[neighbor];

        for (int target = 0; target < 4; ++target) {
            scatter(own, contribution.sites[target], jump * contribution.own_row[target]);
            scatter(neighbor, contribution.sites[target], -jump * contribution.neighbor_row[target]);
        }
    }

    return blocks;
}

// Everything is expressed in the representative cell's chart, where the
// edge, both sites and both opposite sites were reported together, so every
// difference below is chart-consistent.
template<fields::flat::Field FieldType>
typename CapacityHessian<FieldType>::EdgeContribution CapacityHessian<FieldType>::edge_contribution(
    const Diagram& diagram,
    size_t cell,
    const CellEdgeInfo& edge
) const {
    EdgeContribution contribution;
    contribution.sites = {
        cell,
        edge.neighbor_index,
        pinned_site(edge.source_cut, cell),
        pinned_site(edge.target_cut, cell)
    };
    contribution.own_row = {Matrix3::Zero(), Matrix3::Zero(), Matrix3::Zero(), Matrix3::Zero()};
    contribution.neighbor_row = {Matrix3::Zero(), Matrix3::Zero(), Matrix3::Zero(), Matrix3::Zero()};
    contribution.active = false;

    if (cell == edge.neighbor_index) {
        return contribution;
    }

    Vector2 own_site = diagram.site(cell);
    Vector2 neighbor_site = edge.neighbor_position;
    double separation = (neighbor_site - own_site).norm();

    if (separation < GEOMETRIC_EPSILON || edge.boundary.length() < GEOMETRIC_EPSILON) {
        return contribution;
    }

    Vector3 a = planar(own_site);
    Vector3 b = planar(neighbor_site);
    Vector3 normal = (b - a) / separation;

    auto integrals = _field.integrals(edge.boundary);
    double mass = integrals.mass;
    Vector3 moment = integrals.first_moment;
    Vector3 gradient_mass = _field.gradient_masses(edge.boundary);
    Matrix3 gradient_first = _field.gradient_first_moments(edge.boundary);
    std::array<Matrix3, 3> gradient_second = _field.gradient_second_moments(edge.boundary);

    double g0 = normal.dot(gradient_mass);
    Vector3 g1 = gradient_first * normal;
    Matrix3 g2 = normal.x() * gradient_second[0] + normal.y() * gradient_second[1];

    // The sweep: the line's normal displacement is affine along the edge,
    // and a straight line neither stretches nor bends under it, so the
    // integrals respond through the density's normal derivative and the
    // sweep of the integrand itself.
    Response response;
    response.mass = {
        Eigen::RowVector3d((g1 - g0 * a).transpose()) / separation,
        Eigen::RowVector3d(-(g1 - g0 * b).transpose()) / separation,
        Eigen::RowVector3d::Zero(),
        Eigen::RowVector3d::Zero()
    };
    response.moment = {
        (g2 + normal * moment.transpose() - (g1 + normal * mass) * a.transpose()) / separation,
        -(g2 + normal * moment.transpose() - (g1 + normal * mass) * b.transpose()) / separation,
        Matrix3::Zero(),
        Matrix3::Zero()
    };

    Vector3 tangent = planar(edge.boundary.direction()) / edge.boundary.length();

    add_endpoint(response, edge.boundary.source(), own_site, neighbor_site, edge.source_cut, 2, tangent, -1.0);
    add_endpoint(response, edge.boundary.target(), own_site, neighbor_site, edge.target_cut, 3, tangent, 1.0);

    // The sweep rates about each site, F = (moment - site mass) / sep;
    // their derivatives collect the integral responses, the direct site
    // term and the quotient's separation term.
    Vector3 sweep_own = (moment - a * mass) / separation;
    Vector3 sweep_neighbor = (moment - b * mass) / separation;

    for (int target = 0; target < 4; ++target) {
        contribution.own_row[target] = (response.moment[target] - a * response.mass[target]) / separation;
        contribution.neighbor_row[target] = (response.moment[target] - b * response.mass[target]) / separation;
    }

    contribution.own_row[0] += (-mass * planar_identity() + sweep_own * normal.transpose()) / separation;
    contribution.own_row[1] -= sweep_own * normal.transpose() / separation;
    contribution.neighbor_row[1] += (-mass * planar_identity() - sweep_neighbor * normal.transpose()) / separation;
    contribution.neighbor_row[0] += sweep_neighbor * normal.transpose() / separation;

    contribution.active = true;
    return contribution;
}

// The endpoint is equidistant from the two edge sites and pinned by one
// more equation -- equidistance from a third site, or the wall it lies on,
// which does not move -- so its velocity solves the differentiated pair;
// what the integrals feel is its tangential component, an endpoint flux of
// the integrand.
template<fields::flat::Field FieldType>
void CapacityHessian<FieldType>::add_endpoint(
    Response& response,
    const Vector2& endpoint,
    const Vector2& own,
    const Vector2& neighbor,
    const Cut& cut,
    int opposite_slot,
    const Vector3& tangent,
    double sign
) const {
    const Bisector* bisector = std::get_if<Bisector>(&cut);

    Eigen::Matrix2d system;
    system.row(0) = (neighbor - own).transpose();

    if (bisector != nullptr) {
        system.row(1) = (bisector->neighbor_position - own).transpose();
    } else {
        system.row(1) = std::get<Wall>(cut).inward_normal.transpose();
    }

    if (std::abs(system.determinant()) < GEOMETRIC_EPSILON) {
        return;
    }

    Eigen::RowVector2d tangential = Eigen::RowVector2d(tangent.x(), tangent.y()) * system.inverse();

    Eigen::RowVector2d own_row = tangential[0] * (endpoint - own).transpose();
    Eigen::RowVector2d neighbor_row = -tangential[0] * (endpoint - neighbor).transpose();

    std::array<Eigen::RowVector3d, 4> velocity{
        Eigen::RowVector3d::Zero(),
        Eigen::RowVector3d(neighbor_row[0], neighbor_row[1], 0.0),
        Eigen::RowVector3d::Zero(),
        Eigen::RowVector3d::Zero()
    };

    if (bisector != nullptr) {
        own_row += tangential[1] * (endpoint - own).transpose();
        Eigen::RowVector2d opposite_row = -tangential[1] * (endpoint - bisector->neighbor_position).transpose();
        velocity[opposite_slot] = Eigen::RowVector3d(opposite_row[0], opposite_row[1], 0.0);
    }

    velocity[0] = Eigen::RowVector3d(own_row[0], own_row[1], 0.0);

    double density = _field.value(endpoint);
    Vector3 position = planar(endpoint);

    for (int target = 0; target < 4; ++target) {
        response.mass[target] += sign * density * velocity[target];
        response.moment[target] += sign * density * position * velocity[target];
    }
}

template<fields::flat::Field FieldType>
size_t CapacityHessian<FieldType>::pinned_site(const Cut& cut, size_t cell) {
    const Bisector* bisector = std::get_if<Bisector>(&cut);
    return bisector != nullptr ? bisector->neighbor_index : cell;
}

template<fields::flat::Field FieldType>
Matrix3 CapacityHessian<FieldType>::planar_identity() {
    Matrix3 identity = Matrix3::Identity();
    identity(2, 2) = 0.0;
    return identity;
}

template<fields::flat::Field FieldType>
Vector3 CapacityHessian<FieldType>::planar(const Vector2& point) {
    return Vector3(point.x(), point.y(), 0.0);
}

} // namespace geometry_art::voronoi::flat

#endif //GEOMETRY_ART_VORONOI_FLAT_OPTIMIZERS_CAPACITY_HESSIAN_HPP_
