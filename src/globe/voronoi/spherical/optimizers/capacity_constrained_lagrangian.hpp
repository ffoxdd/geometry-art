#ifndef GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_CAPACITY_CONSTRAINED_LAGRANGIAN_HPP_
#define GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_CAPACITY_CONSTRAINED_LAGRANGIAN_HPP_

#include "../../../types.hpp"
#include "../../../fields/spherical/field.hpp"
#include "../../../fields/spherical/polynomial_field.hpp"
#include "../core/sphere.hpp"
#include "../../../fields/region_integrals.hpp"
#include "../../../geometry/spherical/polygon/polygon.hpp"
#include "../../capacity_jacobian.hpp"
#include "../../edge_slots.hpp"
#include "../../state.hpp"
#include "../../../math/polynomial/moments.hpp"
#include "../../../std_ext/parallel_for.hpp"
#include <CGAL/assertions.h>
#include <Eigen/Core>
#include <cmath>
#include <cstddef>
#include <map>
#include <utility>
#include <vector>

namespace globe::voronoi::spherical {

using globe::math::polynomial::Moments;

struct LagrangianEvaluation {
    double value;
    double cvt_energy;
    std::vector<double> capacity_errors;
    std::vector<Vector3> site_gradients;

    [[nodiscard]] double root_mean_square_capacity_error() const;
    [[nodiscard]] double max_absolute_capacity_error() const;
};

template<fields::spherical::Field FieldType = fields::spherical::PolynomialField>
class CapacityConstrainedLagrangian {
 public:
    CapacityConstrainedLagrangian(FieldType field, double target_mass);

    [[nodiscard]] double target_mass() const { return _target_mass; }

    [[nodiscard]] std::vector<CellState> cell_states(const Sphere& sphere) const;
    [[nodiscard]] DiagramState sphere_state(const Sphere& sphere) const;
    [[nodiscard]] LagrangianEvaluation evaluate(
        const Sphere& sphere,
        const std::vector<double>& multipliers,
        double penalty
    ) const;
    [[nodiscard]] LagrangianEvaluation evaluate(
        const Sphere& sphere,
        const DiagramState& state,
        const std::vector<double>& multipliers,
        double penalty
    ) const;

 private:
    struct CellBuild {
        CellState cell;
        std::vector<EdgeState> edges;
    };

    FieldType _field;
    double _target_mass;

    [[nodiscard]] std::vector<Moments> shared_arc_moments(
        const std::vector<std::vector<CellEdgeInfo>>& cell_edges,
        const EdgeSlots& slots
    ) const;
    [[nodiscard]] CellBuild build_cell(
        const std::vector<CellEdgeInfo>& cell_edges,
        const std::vector<size_t>& slots,
        const std::vector<Moments>& shared_moments,
        const std::vector<Vector3>& sites,
        size_t own_index
    ) const;

    [[nodiscard]] std::vector<Vector3> site_gradients(
        const Sphere& sphere,
        const DiagramState& state,
        const std::vector<double>& weights
    ) const;
};

inline double LagrangianEvaluation::root_mean_square_capacity_error() const {
    double sum = 0.0;

    for (double error : capacity_errors) {
        sum += error * error;
    }

    return std::sqrt(sum / static_cast<double>(capacity_errors.size()));
}

inline double LagrangianEvaluation::max_absolute_capacity_error() const {
    double maximum = 0.0;

    for (double error : capacity_errors) {
        maximum = std::max(maximum, std::abs(error));
    }

    return maximum;
}

template<fields::spherical::Field FieldType>
CapacityConstrainedLagrangian<FieldType>::CapacityConstrainedLagrangian(FieldType field, double target_mass) :
    _field(std::move(field)),
    _target_mass(target_mass) {
}

template<fields::spherical::Field FieldType>
std::vector<CellState> CapacityConstrainedLagrangian<FieldType>::cell_states(const Sphere& sphere) const {
    return sphere_state(sphere).cells;
}

template<fields::spherical::Field FieldType>
DiagramState CapacityConstrainedLagrangian<FieldType>::sphere_state(const Sphere& sphere) const {
    size_t count = sphere.size();
    std::vector<std::vector<CellEdgeInfo>> cell_edges(count);

    for (size_t k = 0; k < count; ++k) {
        cell_edges[k] = sphere.cell_edges(k);
    }

    EdgeSlots slots = EdgeSlots::build(cell_edges);
    std::vector<Moments> moments = shared_arc_moments(cell_edges, slots);

    std::vector<Vector3> sites(count);

    for (size_t k = 0; k < count; ++k) {
        sites[k] = to_vector3(sphere.site(k));
    }

    DiagramState state;
    state.cells.resize(count);
    state.edges.resize(count);

    std_ext::parallel_for(count, [&](size_t k) {
        CellBuild built = build_cell(cell_edges[k], slots.slots_by_cell[k], moments, sites, k);
        state.cells[k] = built.cell;
        state.edges[k] = std::move(built.edges);
    });

    return state;
}

template<fields::spherical::Field FieldType>
std::vector<Moments> CapacityConstrainedLagrangian<FieldType>::shared_arc_moments(
    const std::vector<std::vector<CellEdgeInfo>>& cell_edges,
    const EdgeSlots& slots
) const {
    int moment_degree = _field.degree() + 1;
    std::vector<Moments> moments(slots.count(), Moments(moment_degree));

    std_ext::parallel_for(slots.count(), [&](size_t slot) {
        auto [cell, position] = slots.representatives[slot];
        moments[slot] = cell_edges[cell][position].arc.moments(moment_degree);
    });

    return moments;
}

template<fields::spherical::Field FieldType>
typename CapacityConstrainedLagrangian<FieldType>::CellBuild
CapacityConstrainedLagrangian<FieldType>::build_cell(
    const std::vector<CellEdgeInfo>& cell_edges,
    const std::vector<size_t>& slots,
    const std::vector<Moments>& shared_moments,
    const std::vector<Vector3>& sites,
    size_t own_index
) const {
    std::vector<Arc> arcs;
    std::vector<Moments> arc_moments;
    CellBuild built;
    arcs.reserve(cell_edges.size());
    arc_moments.reserve(cell_edges.size());
    built.edges.reserve(cell_edges.size());
    const Vector3& own = sites[own_index];

    for (size_t position = 0; position < cell_edges.size(); ++position) {
        const CellEdgeInfo& edge = cell_edges[position];
        const Moments& moments = shared_moments[slots[position]];
        const Vector3& neighbor = sites[edge.neighbor_index];

        arcs.push_back(edge.arc);
        arc_moments.push_back(moments);
        auto integrals = _field.integrals(edge.arc, moments);
        built.edges.push_back(EdgeState{
            edge.neighbor_index,
            (neighbor - own).norm(),
            integrals.mass,
            integrals.first_moment - own * integrals.mass,
            integrals.first_moment - neighbor * integrals.mass
        });
    }

    // On the unit sphere the squared-norm moment is the mass itself.
    auto integrals = _field.integrals(Polygon(arcs), arc_moments);
    built.cell = CellState{integrals.mass, integrals.first_moment, integrals.mass};

    return built;
}

template<fields::spherical::Field FieldType>
LagrangianEvaluation CapacityConstrainedLagrangian<FieldType>::evaluate(
    const Sphere& sphere,
    const std::vector<double>& multipliers,
    double penalty
) const {
    return evaluate(sphere, sphere_state(sphere), multipliers, penalty);
}

template<fields::spherical::Field FieldType>
LagrangianEvaluation CapacityConstrainedLagrangian<FieldType>::evaluate(
    const Sphere& sphere,
    const DiagramState& state,
    const std::vector<double>& multipliers,
    double penalty
) const {
    size_t count = sphere.size();
    CGAL_precondition(multipliers.size() == count);

    const std::vector<CellState>& states = state.cells;
    std::vector<double> capacity_errors(count);
    std::vector<double> weights(count);
    double cvt_energy = 0.0;
    double value = 0.0;

    for (size_t i = 0; i < count; ++i) {
        Vector3 site = to_vector3(sphere.site(i));
        capacity_errors[i] = states[i].mass - _target_mass;
        weights[i] = multipliers[i] + penalty * capacity_errors[i];
        cvt_energy += states[i].squared_norm_moment -
            2.0 * site.dot(states[i].first_moment) + site.squaredNorm() * states[i].mass;
        value += multipliers[i] * capacity_errors[i] + 0.5 * penalty * capacity_errors[i] * capacity_errors[i];
    }

    value += cvt_energy;

    return LagrangianEvaluation{
        value,
        cvt_energy,
        std::move(capacity_errors),
        site_gradients(sphere, state, weights)
    };
}

template<fields::spherical::Field FieldType>
std::vector<Vector3> CapacityConstrainedLagrangian<FieldType>::site_gradients(
    const Sphere& sphere,
    const DiagramState& state,
    const std::vector<double>& weights
) const {
    size_t count = sphere.size();
    std::vector<Vector3> gradients = CapacityJacobian(state).transpose_apply(weights);

    for (size_t k = 0; k < count; ++k) {
        gradients[k] -= 2.0 * state.cells[k].first_moment;
    }

    return gradients;
}

} // namespace globe::voronoi::spherical

#endif //GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_CAPACITY_CONSTRAINED_LAGRANGIAN_HPP_
