#ifndef GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_CAPACITY_CONSTRAINED_LAGRANGIAN_HPP_
#define GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_CAPACITY_CONSTRAINED_LAGRANGIAN_HPP_

#include "../../../types.hpp"
#include "../../../fields/spherical/field.hpp"
#include "../../../fields/spherical/polynomial_field.hpp"
#include "../core/sphere.hpp"
#include "../../../fields/spherical/region_integrals.hpp"
#include "../../../geometry/spherical/polygon/polygon.hpp"
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

struct CellState {
    double mass;
    Vector3 first_moment;
};

struct EdgeState {
    size_t neighbor_index;
    fields::spherical::RegionIntegrals integrals;
};

struct SphereState {
    std::vector<CellState> cells;
    std::vector<std::vector<EdgeState>> edges;
};

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
    [[nodiscard]] SphereState sphere_state(const Sphere& sphere) const;
    [[nodiscard]] LagrangianEvaluation evaluate(
        const Sphere& sphere,
        const std::vector<double>& multipliers,
        double penalty
    ) const;

 private:
    using EdgeKey = std::pair<size_t, size_t>;

    struct SharedArcs {
        std::vector<Moments> moments;
        std::vector<std::vector<size_t>> slots_by_cell;
    };

    struct CellBuild {
        CellState cell;
        std::vector<EdgeState> edges;
    };

    FieldType _field;
    double _target_mass;

    [[nodiscard]] SharedArcs shared_arcs(const std::vector<std::vector<CellEdgeInfo>>& cell_edges) const;
    [[nodiscard]] CellBuild build_cell(
        const std::vector<CellEdgeInfo>& cell_edges,
        const std::vector<size_t>& slots,
        const std::vector<Moments>& shared_moments
    ) const;

    [[nodiscard]] std::vector<Vector3> site_gradients(
        const Sphere& sphere,
        const SphereState& state,
        const std::vector<double>& weights
    ) const;
    [[nodiscard]] static Vector3 edge_mass_gradient(
        const Vector3& site,
        const Vector3& neighbor,
        const fields::spherical::RegionIntegrals& edge_integrals
    );
    [[nodiscard]] static EdgeKey edge_key(size_t a, size_t b);
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
SphereState CapacityConstrainedLagrangian<FieldType>::sphere_state(const Sphere& sphere) const {
    size_t count = sphere.size();
    std::vector<std::vector<CellEdgeInfo>> cell_edges(count);

    for (size_t k = 0; k < count; ++k) {
        cell_edges[k] = sphere.cell_edges(k);
    }

    SharedArcs shared = shared_arcs(cell_edges);
    SphereState state;
    state.cells.resize(count);
    state.edges.resize(count);

    std_ext::parallel_for(count, [&](size_t k) {
        CellBuild built = build_cell(cell_edges[k], shared.slots_by_cell[k], shared.moments);
        state.cells[k] = built.cell;
        state.edges[k] = std::move(built.edges);
    });

    return state;
}

template<fields::spherical::Field FieldType>
typename CapacityConstrainedLagrangian<FieldType>::SharedArcs
CapacityConstrainedLagrangian<FieldType>::shared_arcs(
    const std::vector<std::vector<CellEdgeInfo>>& cell_edges
) const {
    std::map<EdgeKey, size_t> slot_by_key;
    std::vector<Arc> arcs;
    SharedArcs shared;
    shared.slots_by_cell.resize(cell_edges.size());

    for (size_t k = 0; k < cell_edges.size(); ++k) {
        shared.slots_by_cell[k].reserve(cell_edges[k].size());

        for (const CellEdgeInfo& edge : cell_edges[k]) {
            auto [iterator, inserted] = slot_by_key.try_emplace(edge_key(k, edge.neighbor_index), arcs.size());
            if (inserted) {
                arcs.push_back(edge.arc);
            }

            shared.slots_by_cell[k].push_back(iterator->second);
        }
    }

    int moment_degree = _field.degree() + 1;
    shared.moments.assign(arcs.size(), Moments(moment_degree));

    std_ext::parallel_for(arcs.size(), [&](size_t slot) {
        shared.moments[slot] = arcs[slot].moments(moment_degree);
    });

    return shared;
}

template<fields::spherical::Field FieldType>
typename CapacityConstrainedLagrangian<FieldType>::CellBuild
CapacityConstrainedLagrangian<FieldType>::build_cell(
    const std::vector<CellEdgeInfo>& cell_edges,
    const std::vector<size_t>& slots,
    const std::vector<Moments>& shared_moments
) const {
    std::vector<Arc> arcs;
    std::vector<Moments> arc_moments;
    CellBuild built;
    arcs.reserve(cell_edges.size());
    arc_moments.reserve(cell_edges.size());
    built.edges.reserve(cell_edges.size());

    for (size_t position = 0; position < cell_edges.size(); ++position) {
        const CellEdgeInfo& edge = cell_edges[position];
        const Moments& moments = shared_moments[slots[position]];

        arcs.push_back(edge.arc);
        arc_moments.push_back(moments);
        built.edges.push_back(EdgeState{edge.neighbor_index, _field.integrals(edge.arc, moments)});
    }

    auto integrals = _field.integrals(Polygon(arcs), arc_moments);
    built.cell = CellState{integrals.mass, integrals.first_moment};

    return built;
}

template<fields::spherical::Field FieldType>
LagrangianEvaluation CapacityConstrainedLagrangian<FieldType>::evaluate(
    const Sphere& sphere,
    const std::vector<double>& multipliers,
    double penalty
) const {
    size_t count = sphere.size();
    CGAL_precondition(multipliers.size() == count);

    SphereState state = sphere_state(sphere);
    const std::vector<CellState>& states = state.cells;
    std::vector<double> capacity_errors(count);
    std::vector<double> weights(count);
    double cvt_energy = 0.0;
    double value = 0.0;

    for (size_t i = 0; i < count; ++i) {
        Vector3 site = to_vector3(sphere.site(i));
        capacity_errors[i] = states[i].mass - _target_mass;
        weights[i] = multipliers[i] + penalty * capacity_errors[i];
        cvt_energy += 2.0 * (states[i].mass - site.dot(states[i].first_moment));
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
    const SphereState& state,
    const std::vector<double>& weights
) const {
    size_t count = sphere.size();
    std::vector<Vector3> sites(count);
    for (size_t k = 0; k < count; ++k) {
        sites[k] = to_vector3(sphere.site(k));
    }

    std::vector<Vector3> gradients(count);

    for (size_t k = 0; k < count; ++k) {
        Vector3 gradient = -2.0 * state.cells[k].first_moment;

        for (const EdgeState& edge : state.edges[k]) {
            gradient += (weights[k] - weights[edge.neighbor_index]) *
                edge_mass_gradient(sites[k], sites[edge.neighbor_index], edge.integrals);
        }

        gradients[k] = gradient;
    }

    return gradients;
}

template<fields::spherical::Field FieldType>
Vector3 CapacityConstrainedLagrangian<FieldType>::edge_mass_gradient(
    const Vector3& site,
    const Vector3& neighbor,
    const fields::spherical::RegionIntegrals& edge_integrals
) {
    double separation = (neighbor - site).norm();

    if (separation < GEOMETRIC_EPSILON) {
        return Vector3::Zero();
    }

    return (edge_integrals.first_moment - site * edge_integrals.mass) / separation;
}

template<fields::spherical::Field FieldType>
typename CapacityConstrainedLagrangian<FieldType>::EdgeKey
CapacityConstrainedLagrangian<FieldType>::edge_key(size_t a, size_t b) {
    return a < b ? EdgeKey{a, b} : EdgeKey{b, a};
}

} // namespace globe::voronoi::spherical

#endif //GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_CAPACITY_CONSTRAINED_LAGRANGIAN_HPP_
