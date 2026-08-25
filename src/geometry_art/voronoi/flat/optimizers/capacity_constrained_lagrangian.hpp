#ifndef GEOMETRY_ART_VORONOI_FLAT_OPTIMIZERS_CAPACITY_CONSTRAINED_LAGRANGIAN_HPP_
#define GEOMETRY_ART_VORONOI_FLAT_OPTIMIZERS_CAPACITY_CONSTRAINED_LAGRANGIAN_HPP_

#include "../core/periodic_slots.hpp"
#include "../core/torus.hpp"
#include "../../capacity_jacobian.hpp"
#include "../../lagrangian_evaluation.hpp"
#include "../../state.hpp"
#include "../../../fields/flat/field.hpp"
#include "../../../std_ext/parallel_for.hpp"
#include "../../../types.hpp"
#include <CGAL/assertions.h>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace geometry_art::voronoi::flat {

// The augmented Lagrangian of capacity-constrained CVT on the flat torus.
//
// Bisectors are shared through slots keyed by the site pair AND the period
// offset between the two charts: a torus can hold two distinct bisectors
// between the same two sites, one across each seam, and a cell can even
// border itself. Everything stored is chart-invariant, so both cells read
// one computation.
template<fields::flat::Field FieldType>
class CapacityConstrainedLagrangian {
 public:
    CapacityConstrainedLagrangian(FieldType field, double target_mass);

    [[nodiscard]] double target_mass() const { return _target_mass; }

    [[nodiscard]] DiagramState diagram_state(const Torus& torus) const;
    [[nodiscard]] LagrangianEvaluation evaluate(
        const Torus& torus,
        const DiagramState& state,
        const std::vector<double>& multipliers,
        double penalty
    ) const;
    [[nodiscard]] LagrangianEvaluation evaluate(
        const Torus& torus,
        const std::vector<double>& multipliers,
        double penalty
    ) const;

 private:
    struct SlotState {
        size_t source_cell;
        double separation;
        double mass;
        Vector3 moment_about_source;
        Vector3 moment_about_target;
    };

    FieldType _field;
    double _target_mass;

    [[nodiscard]] static Vector3 planar(const Vector2& point);
};

template<fields::flat::Field FieldType>
CapacityConstrainedLagrangian<FieldType>::CapacityConstrainedLagrangian(FieldType field, double target_mass) :
    _field(std::move(field)),
    _target_mass(target_mass) {
}

template<fields::flat::Field FieldType>
DiagramState CapacityConstrainedLagrangian<FieldType>::diagram_state(const Torus& torus) const {
    size_t count = torus.size();
    std::vector<std::vector<CellEdgeInfo>> cell_edges(count);

    for (size_t k = 0; k < count; ++k) {
        cell_edges[k] = torus.cell_edges(k);
    }

    PeriodicSlots shared = PeriodicSlots::build(torus, cell_edges);
    std::vector<SlotState> slots(shared.count());

    std_ext::parallel_for(shared.count(), [&](size_t slot) {
        auto [cell, position] = shared.representatives[slot];
        const CellEdgeInfo& edge = cell_edges[cell][position];
        auto integrals = _field.integrals(edge.boundary);
        Vector3 own = planar(torus.site(cell));
        Vector3 neighbor = planar(edge.neighbor_position);

        slots[slot] = SlotState{
            cell,
            (neighbor - own).norm(),
            integrals.mass,
            integrals.first_moment - own * integrals.mass,
            integrals.first_moment - neighbor * integrals.mass
        };
    });

    DiagramState state;
    state.cells.resize(count);
    state.edges.resize(count);

    std_ext::parallel_for(count, [&](size_t k) {
        Polygon cell = torus.cell(k);
        auto integrals = _field.integrals(cell);
        state.cells[k] = CellState{
            integrals.mass,
            integrals.first_moment,
            _field.squared_norm_moment(cell)
        };

        state.edges[k].reserve(cell_edges[k].size());

        for (size_t position = 0; position < cell_edges[k].size(); ++position) {
            const SlotState& slot = slots[shared.slots_by_cell[k][position]];
            bool from_source = slot.source_cell == k;

            state.edges[k].push_back(EdgeState{
                cell_edges[k][position].neighbor_index,
                slot.separation,
                slot.mass,
                from_source ? slot.moment_about_source : slot.moment_about_target,
                from_source ? slot.moment_about_target : slot.moment_about_source
            });
        }
    });

    return state;
}

template<fields::flat::Field FieldType>
LagrangianEvaluation CapacityConstrainedLagrangian<FieldType>::evaluate(
    const Torus& torus,
    const std::vector<double>& multipliers,
    double penalty
) const {
    return evaluate(torus, diagram_state(torus), multipliers, penalty);
}

template<fields::flat::Field FieldType>
LagrangianEvaluation CapacityConstrainedLagrangian<FieldType>::evaluate(
    const Torus& torus,
    const DiagramState& state,
    const std::vector<double>& multipliers,
    double penalty
) const {
    size_t count = torus.size();
    CGAL_precondition(multipliers.size() == count);

    std::vector<double> capacity_errors(count);
    std::vector<double> weights(count);
    double cvt_energy = 0.0;
    double value = 0.0;

    for (size_t i = 0; i < count; ++i) {
        Vector3 site = torus.site_vector(i);
        capacity_errors[i] = state.cells[i].mass - _target_mass;
        weights[i] = multipliers[i] + penalty * capacity_errors[i];
        cvt_energy += state.cells[i].squared_norm_moment -
            2.0 * site.dot(state.cells[i].first_moment) + site.squaredNorm() * state.cells[i].mass;
        value += multipliers[i] * capacity_errors[i] + 0.5 * penalty * capacity_errors[i] * capacity_errors[i];
    }

    value += cvt_energy;

    std::vector<Vector3> gradients = CapacityJacobian(state).transpose_apply(weights);

    // The CVT gradient on a flat domain keeps its position term; there is
    // no radial direction for a tangent projection to discard it into.
    for (size_t k = 0; k < count; ++k) {
        gradients[k] += 2.0 * (state.cells[k].mass * torus.site_vector(k) - state.cells[k].first_moment);
    }

    return LagrangianEvaluation{value, cvt_energy, std::move(capacity_errors), std::move(gradients)};
}

template<fields::flat::Field FieldType>
Vector3 CapacityConstrainedLagrangian<FieldType>::planar(const Vector2& point) {
    return Vector3(point.x(), point.y(), 0.0);
}

} // namespace geometry_art::voronoi::flat

#endif //GEOMETRY_ART_VORONOI_FLAT_OPTIMIZERS_CAPACITY_CONSTRAINED_LAGRANGIAN_HPP_
