#ifndef GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_CAPACITY_HESSIAN_HPP_
#define GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_CAPACITY_HESSIAN_HPP_

#include "cvt_hessian.hpp"
#include "../../state.hpp"
#include "../../../types.hpp"
#include "../../../fields/spherical/field.hpp"
#include "../../../fields/spherical/polynomial_field.hpp"
#include "../../../std_ext/parallel_for.hpp"
#include "../core/sphere.hpp"
#include <CGAL/assertions.h>
#include <Eigen/Dense>
#include <array>
#include <cstddef>
#include <utility>
#include <vector>

namespace globe::voronoi::spherical {

// Second derivatives of the weighted cell masses with respect to the sites:
// the curvature of `sum_i w_i c_i`, which is the term the Gauss-Newton model
// of the augmented Lagrangian leaves out. It does not vanish at a solution,
// because mass has a price wherever the density varies.
//
// A weighted mass changes only where bisectors sweep, so the curvature is a
// sum over bisector arcs, and differentiating one arc's sweep splits into
// three motions:
//
//   - the great circle rotates as its two sites move, felt through the
//     density and its gradient along the arc -- arc moments the field
//     already knows how to integrate exactly;
//   - the separation of the sites scales the sweep rate;
//   - the endpoints slide, each the meeting point of three bisectors, with
//     velocities given by the implicit function theorem on the three
//     equidistance equations.
//
// The endpoint of an arc is shared with one more site, so blocks couple a
// site to every Delaunay neighbour and to nothing else: the sparsity is the
// energy Hessian's.
template<fields::spherical::Field FieldType = fields::spherical::PolynomialField>
class CapacityHessian {
 public:
    explicit CapacityHessian(FieldType field);

    [[nodiscard]] HessianBlocks assemble(
        const Sphere& sphere,
        const DiagramState& state,
        const std::vector<Vector3>& sites,
        const std::vector<double>& weights
    ) const;

 private:
    struct EdgeContribution {
        std::array<size_t, 4> sites;
        std::array<Matrix3, 4> blocks;
        double gauge;
    };

    FieldType _field;

    [[nodiscard]] EdgeContribution edge_contribution(
        const CellEdgeInfo& edge,
        const fields::RegionIntegrals& integrals,
        size_t own_index,
        const std::vector<Vector3>& sites
    ) const;

    static void add_endpoint(
        EdgeContribution& contribution,
        int opposite_slot,
        const Vector3& vertex,
        const Vector3& outward_tangent,
        double density,
        const std::vector<Vector3>& sites,
        const Matrix3& rotation_rows,
        double separation
    );

    [[nodiscard]] static Matrix3 rotation_curvature(
        const std::array<Matrix3, 3>& gradient_second_moments,
        const Vector3& first_moment
    );

    [[nodiscard]] static Matrix3 cross_matrix(const Vector3& vector);
};

template<fields::spherical::Field FieldType>
CapacityHessian<FieldType>::CapacityHessian(FieldType field) :
    _field(std::move(field)) {
}

template<fields::spherical::Field FieldType>
HessianBlocks CapacityHessian<FieldType>::assemble(
    const Sphere& sphere,
    const DiagramState& state,
    const std::vector<Vector3>& sites,
    const std::vector<double>& weights
) const {
    size_t count = sphere.size();
    std::vector<std::vector<CellEdgeInfo>> cell_edges(count);

    for (size_t k = 0; k < count; ++k) {
        cell_edges[k] = sphere.cell_edges(k);
    }

    std::vector<std::pair<size_t, size_t>> unique_edges;

    for (size_t k = 0; k < count; ++k) {
        for (size_t position = 0; position < cell_edges[k].size(); ++position) {
            if (cell_edges[k][position].neighbor_index > k) {
                unique_edges.emplace_back(k, position);
            }
        }
    }

    std::vector<EdgeContribution> contributions(unique_edges.size());

    std_ext::parallel_for(unique_edges.size(), [&](size_t slot) {
        auto [k, position] = unique_edges[slot];
        CGAL_precondition(state.edges[k][position].neighbor_index == cell_edges[k][position].neighbor_index);

        contributions[slot] = edge_contribution(
            cell_edges[k][position],
            state.edges[k][position].integrals,
            k,
            sites
        );
    });

    HessianBlocks result;
    result.diagonal.assign(count, Matrix3::Zero());
    result.neighbors.resize(count);

    for (size_t k = 0; k < count; ++k) {
        result.neighbors[k].reserve(cell_edges[k].size());

        for (const CellEdgeInfo& edge : cell_edges[k]) {
            result.neighbors[k].push_back(NeighborBlock{edge.neighbor_index, Matrix3::Zero()});
        }
    }

    auto block_of = [&](size_t row, size_t column) -> Matrix3& {
        if (row == column) {
            return result.diagonal[row];
        }

        for (NeighborBlock& block : result.neighbors[row]) {
            if (block.neighbor_index == column) {
                return block.value;
            }
        }

        CGAL_error_msg("a capacity curvature block reached beyond the Delaunay neighbours");
        return result.diagonal[row];
    };

    for (size_t slot = 0; slot < unique_edges.size(); ++slot) {
        auto [k, position] = unique_edges[slot];
        size_t j = cell_edges[k][position].neighbor_index;
        const EdgeContribution& contribution = contributions[slot];
        double jump = weights[k] - weights[j];

        for (int entry = 0; entry < 4; ++entry) {
            size_t column = contribution.sites[entry];
            block_of(k, column) += jump * contribution.blocks[entry];
            block_of(j, column) -= jump * contribution.blocks[entry];
        }

        // The stored sweep rates are expressed about each site, and the
        // curvature of that expression is what keeps this assembly
        // consistent with the gradient the optimizer actually descends.
        result.diagonal[k] -= jump * contribution.gauge * Matrix3::Identity();
        result.diagonal[j] += jump * contribution.gauge * Matrix3::Identity();
    }

    return result;
}

// The derivative of one arc's sweep rate `F = (integral of density times
// position) / separation` along every site motion it depends on, before the
// weight jump is applied.
template<fields::spherical::Field FieldType>
typename CapacityHessian<FieldType>::EdgeContribution CapacityHessian<FieldType>::edge_contribution(
    const CellEdgeInfo& edge,
    const fields::RegionIntegrals& integrals,
    size_t own_index,
    const std::vector<Vector3>& sites
) const {
    EdgeContribution contribution{
        {own_index, edge.neighbor_index, edge.source_opposite_index, edge.target_opposite_index},
        {Matrix3::Zero(), Matrix3::Zero(), Matrix3::Zero(), Matrix3::Zero()},
        0.0
    };

    const Vector3& own = sites[own_index];
    const Vector3& neighbor = sites[edge.neighbor_index];
    Vector3 difference = own - neighbor;
    double separation = difference.norm();

    if (separation < GEOMETRIC_EPSILON) {
        return contribution;
    }

    Vector3 normal = difference / separation;
    Matrix3 normal_cross = cross_matrix(normal);

    // The bisector rotates with the minimal angular velocity taking its
    // normal along, so a site motion turns into `omega = normal x (u_own -
    // u_neighbor) / separation`, and these rows convert curvature in omega
    // into curvature in the site motions.
    Matrix3 rotation_rows = normal_cross / separation;

    Matrix3 swept = rotation_curvature(_field.gradient_second_moments(edge.arc), integrals.first_moment) *
        rotation_rows / separation;

    Matrix3 metric = integrals.first_moment * normal.transpose() / (separation * separation);

    contribution.blocks[0] = swept - metric;
    contribution.blocks[1] = -contribution.blocks[0];
    contribution.gauge = integrals.mass / separation;

    add_endpoint(
        contribution,
        2,
        edge.arc.source(),
        -edge.arc.normal().cross(edge.arc.source()),
        _field.value(edge.arc.source()),
        sites,
        rotation_rows,
        separation
    );

    add_endpoint(
        contribution,
        3,
        edge.arc.target(),
        edge.arc.normal().cross(edge.arc.target()),
        _field.value(edge.arc.target()),
        sites,
        rotation_rows,
        separation
    );

    return contribution;
}

// An endpoint is a Voronoi vertex: equidistant from the edge's two sites and
// one more. Its velocity solves the differentiated equidistance equations,
// and the arc gains or loses length at the rate the endpoint slides along
// the rotating circle.
template<fields::spherical::Field FieldType>
void CapacityHessian<FieldType>::add_endpoint(
    EdgeContribution& contribution,
    int opposite_slot,
    const Vector3& vertex,
    const Vector3& outward_tangent,
    double density,
    const std::vector<Vector3>& sites,
    const Matrix3& rotation_rows,
    double separation
) {
    const Vector3& own = sites[contribution.sites[0]];
    const Vector3& neighbor = sites[contribution.sites[1]];
    const Vector3& opposite = sites[contribution.sites[opposite_slot]];

    Matrix3 equidistance;
    equidistance.row(0) = (own - neighbor).transpose();
    equidistance.row(1) = (own - opposite).transpose();
    equidistance.row(2) = vertex.transpose();

    Eigen::FullPivLU<Matrix3> decomposition(equidistance);

    if (!decomposition.isInvertible()) {
        return;
    }

    Matrix3 inverse = decomposition.inverse();
    Vector3 own_response = inverse * Vector3(-1.0, -1.0, 0.0);
    Vector3 neighbor_response = inverse * Vector3(1.0, 0.0, 0.0);
    Vector3 opposite_response = inverse * Vector3(0.0, 1.0, 0.0);

    Eigen::RowVector3d relative_rows = outward_tangent.transpose() * cross_matrix(vertex) * rotation_rows;

    double factor = density / separation;
    Vector3 scaled_vertex = factor * vertex;

    contribution.blocks[0] +=
        scaled_vertex * (outward_tangent.dot(own_response) * vertex.transpose() + relative_rows);
    contribution.blocks[1] +=
        scaled_vertex * (outward_tangent.dot(neighbor_response) * vertex.transpose() - relative_rows);
    contribution.blocks[opposite_slot] +=
        scaled_vertex * (outward_tangent.dot(opposite_response) * vertex.transpose());
}

// The interior of the arc turns rigidly with the circle, so the integrand's
// change is the density gradient read along the rotation plus the rotation
// of the position factor itself.
template<fields::spherical::Field FieldType>
Matrix3 CapacityHessian<FieldType>::rotation_curvature(
    const std::array<Matrix3, 3>& gradient_second_moments,
    const Vector3& first_moment
) {
    Matrix3 gradient_term;

    for (int axis = 0; axis < 3; ++axis) {
        gradient_term.col(axis) =
            gradient_second_moments[(axis + 2) % 3].col((axis + 1) % 3) -
            gradient_second_moments[(axis + 1) % 3].col((axis + 2) % 3);
    }

    return gradient_term - cross_matrix(first_moment);
}

template<fields::spherical::Field FieldType>
Matrix3 CapacityHessian<FieldType>::cross_matrix(const Vector3& vector) {
    Matrix3 result;

    result <<
        0.0, -vector.z(), vector.y(),
        vector.z(), 0.0, -vector.x(),
        -vector.y(), vector.x(), 0.0;

    return result;
}

} // namespace globe::voronoi::spherical

#endif //GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_CAPACITY_HESSIAN_HPP_
