#ifndef GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_CVT_HESSIAN_HPP_
#define GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_CVT_HESSIAN_HPP_

#include "../../../types.hpp"
#include "../../edge_slots.hpp"
#include "../../hessian_blocks.hpp"
#include "../../../fields/spherical/field.hpp"
#include "../../../fields/spherical/polynomial_field.hpp"
#include "../../../math/normalization.hpp"
#include "../../../std_ext/parallel_for.hpp"
#include "../core/sphere.hpp"
#include <cstddef>
#include <map>
#include <utility>
#include <vector>

namespace globe::voronoi::spherical {

// Second derivatives of the CVT energy with respect to the sites as ambient
// vectors, one dense block per cell plus one per shared bisector.
template<fields::spherical::Field FieldType = fields::spherical::PolynomialField>
class CvtHessian {
 public:
    explicit CvtHessian(FieldType field);

    [[nodiscard]] HessianBlocks assemble(const Sphere& sphere) const;

 private:
    FieldType _field;

    [[nodiscard]] std::vector<Matrix3> shared_second_moments(
        const std::vector<std::vector<CellEdgeInfo>>& cell_edges,
        const EdgeSlots& slots
    ) const;
};

template<fields::spherical::Field FieldType>
CvtHessian<FieldType>::CvtHessian(FieldType field) :
    _field(std::move(field)) {
}

template<fields::spherical::Field FieldType>
HessianBlocks CvtHessian<FieldType>::assemble(const Sphere& sphere) const {
    size_t count = sphere.size();
    std::vector<Vector3> sites(count);
    std::vector<std::vector<CellEdgeInfo>> cell_edges(count);

    for (size_t k = 0; k < count; ++k) {
        sites[k] = to_vector3(sphere.site(k));
        cell_edges[k] = sphere.cell_edges(k);
    }

    EdgeSlots slots = EdgeSlots::build(cell_edges);
    std::vector<Matrix3> second_moments = shared_second_moments(cell_edges, slots);

    HessianBlocks blocks;
    blocks.diagonal.assign(count, Matrix3::Zero());
    blocks.neighbors.resize(count);

    std_ext::parallel_for(count, [&](size_t k) {
        blocks.neighbors[k].reserve(cell_edges[k].size());

        for (size_t position = 0; position < cell_edges[k].size(); ++position) {
            const CellEdgeInfo& edge = cell_edges[k][position];
            double separation = (sites[edge.neighbor_index] - sites[k]).norm();

            if (separation < GEOMETRIC_EPSILON) {
                blocks.neighbors[k].push_back(NeighborBlock{edge.neighbor_index, Matrix3::Zero()});
                continue;
            }

            Matrix3 block = 2.0 * second_moments[slots.slots_by_cell[k][position]] / separation;
            blocks.diagonal[k] -= block;
            blocks.neighbors[k].push_back(NeighborBlock{edge.neighbor_index, block});
        }
    });

    return blocks;
}

// Moving a site sweeps the shared bisector at a rate proportional to the
// point's projection on the displacement, so every block needs the density's
// second moment along one bisector -- the same integral from either side.
template<fields::spherical::Field FieldType>
std::vector<Matrix3> CvtHessian<FieldType>::shared_second_moments(
    const std::vector<std::vector<CellEdgeInfo>>& cell_edges,
    const EdgeSlots& slots
) const {
    std::vector<Matrix3> second_moments(slots.count(), Matrix3::Zero());

    std_ext::parallel_for(slots.count(), [&](size_t slot) {
        auto [cell, position] = slots.representatives[slot];
        second_moments[slot] = _field.second_moment(cell_edges[cell][position].arc);
    });

    return second_moments;
}

} // namespace globe::voronoi::spherical

#endif //GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_CVT_HESSIAN_HPP_
