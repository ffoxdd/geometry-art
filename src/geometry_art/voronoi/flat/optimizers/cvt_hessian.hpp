#ifndef GEOMETRY_ART_VORONOI_FLAT_OPTIMIZERS_CVT_HESSIAN_HPP_
#define GEOMETRY_ART_VORONOI_FLAT_OPTIMIZERS_CVT_HESSIAN_HPP_

#include "../core/periodic_slots.hpp"
#include "../core/torus.hpp"
#include "../../hessian_blocks.hpp"
#include "../../../fields/flat/field.hpp"
#include "../../../std_ext/parallel_for.hpp"
#include "../../../types.hpp"
#include <cstddef>
#include <utility>
#include <vector>

namespace geometry_art::voronoi::flat {

// Second derivatives of the CVT energy on the flat torus. The gradient is
// 2(m s - M1) per cell; differentiating it gives 2mP on the diagonal plus
// the sweep of every bisector, which contributes the relative second
// moments R_ab = integral of rho (x - a)(x - b)^T along the edge -- all
// chart-invariant, so one computation per physical bisector serves both
// sides.
template<fields::flat::Field FieldType>
class CvtHessian {
 public:
    explicit CvtHessian(FieldType field);

    [[nodiscard]] HessianBlocks assemble(const Torus& torus) const;

 private:
    struct SlotCurvature {
        size_t source_cell;
        double separation;
        Matrix3 about_source;
        Matrix3 about_target;
        Matrix3 source_target;
    };

    FieldType _field;

    [[nodiscard]] static Matrix3 relative_second_moment(
        const Matrix3& second_moment,
        const Vector3& first_moment,
        double mass,
        const Vector3& left,
        const Vector3& right
    );
    [[nodiscard]] static Matrix3 planar_identity();
};

template<fields::flat::Field FieldType>
CvtHessian<FieldType>::CvtHessian(FieldType field) :
    _field(std::move(field)) {
}

template<fields::flat::Field FieldType>
HessianBlocks CvtHessian<FieldType>::assemble(const Torus& torus) const {
    size_t count = torus.size();
    std::vector<std::vector<CellEdgeInfo>> cell_edges(count);

    for (size_t k = 0; k < count; ++k) {
        cell_edges[k] = torus.cell_edges(k);
    }

    PeriodicSlots slots = PeriodicSlots::build(torus, cell_edges);
    std::vector<SlotCurvature> curvatures(slots.count());

    std_ext::parallel_for(slots.count(), [&](size_t slot) {
        auto [cell, position] = slots.representatives[slot];
        const CellEdgeInfo& edge = cell_edges[cell][position];
        auto integrals = _field.integrals(edge.boundary);
        Matrix3 second_moment = _field.second_moment(edge.boundary);
        Vector2 own_site = torus.site(cell);
        Vector3 own(own_site.x(), own_site.y(), 0.0);
        Vector3 neighbor(edge.neighbor_position.x(), edge.neighbor_position.y(), 0.0);

        curvatures[slot] = SlotCurvature{
            cell,
            (neighbor - own).norm(),
            relative_second_moment(second_moment, integrals.first_moment, integrals.mass, own, own),
            relative_second_moment(second_moment, integrals.first_moment, integrals.mass, neighbor, neighbor),
            relative_second_moment(second_moment, integrals.first_moment, integrals.mass, own, neighbor)
        };
    });

    HessianBlocks blocks;
    blocks.diagonal.assign(count, Matrix3::Zero());
    blocks.neighbors.resize(count);

    std_ext::parallel_for(count, [&](size_t k) {
        blocks.diagonal[k] = 2.0 * _field.integrals(torus.cell(k)).mass * planar_identity();
        blocks.neighbors[k].reserve(cell_edges[k].size());

        for (size_t position = 0; position < cell_edges[k].size(); ++position) {
            const SlotCurvature& slot = curvatures[slots.slots_by_cell[k][position]];

            if (slot.separation < GEOMETRIC_EPSILON) {
                blocks.neighbors[k].push_back(
                    NeighborBlock{cell_edges[k][position].neighbor_index, Matrix3::Zero()}
                );
                continue;
            }

            bool from_source = slot.source_cell == k;
            const Matrix3& own_own = from_source ? slot.about_source : slot.about_target;
            Matrix3 own_neighbor = from_source ? slot.source_target : slot.source_target.transpose();

            blocks.diagonal[k] -= 2.0 * own_own / slot.separation;
            blocks.neighbors[k].push_back(NeighborBlock{
                cell_edges[k][position].neighbor_index,
                2.0 * own_neighbor / slot.separation
            });
        }
    });

    return blocks;
}

template<fields::flat::Field FieldType>
Matrix3 CvtHessian<FieldType>::relative_second_moment(
    const Matrix3& second_moment,
    const Vector3& first_moment,
    double mass,
    const Vector3& left,
    const Vector3& right
) {
    return second_moment - left * first_moment.transpose() -
        first_moment * right.transpose() + mass * left * right.transpose();
}

template<fields::flat::Field FieldType>
Matrix3 CvtHessian<FieldType>::planar_identity() {
    Matrix3 identity = Matrix3::Identity();
    identity(2, 2) = 0.0;
    return identity;
}

} // namespace geometry_art::voronoi::flat

#endif //GEOMETRY_ART_VORONOI_FLAT_OPTIMIZERS_CVT_HESSIAN_HPP_
