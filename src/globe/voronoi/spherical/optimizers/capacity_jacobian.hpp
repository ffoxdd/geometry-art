#ifndef GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_CAPACITY_JACOBIAN_HPP_
#define GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_CAPACITY_JACOBIAN_HPP_

#include "sphere_state.hpp"
#include "../../../types.hpp"
#include "../../../fields/spherical/region_integrals.hpp"
#include <cstddef>
#include <utility>
#include <vector>

namespace globe::voronoi::spherical {

// Derivatives of every cell's mass with respect to every site. Moving a site
// sweeps only the bisectors it borders, so a cell's mass responds to its own
// site and to its Delaunay neighbours' and to nothing else.
//
// The two directions are kept as operators rather than a matrix: the forward
// one reports how fast each mass changes along a displacement of all sites,
// and the transpose one accumulates a weighted mass gradient, which is
// exactly the constraint part of the Lagrangian's site gradient.
class CapacityJacobian {
 public:
    CapacityJacobian(const SphereState& state, std::vector<Vector3> sites);

    [[nodiscard]] std::vector<double> apply(const std::vector<Vector3>& directions) const;
    [[nodiscard]] std::vector<Vector3> transpose_apply(const std::vector<double>& values) const;

 private:
    struct EdgeGradient {
        size_t neighbor_index;
        Vector3 own;
        Vector3 neighbor;
    };

    std::vector<std::vector<EdgeGradient>> _edges;

    [[nodiscard]] static Vector3 sweep_rate(
        const Vector3& site,
        const Vector3& neighbor,
        const fields::spherical::RegionIntegrals& integrals
    );
};

inline CapacityJacobian::CapacityJacobian(const SphereState& state, std::vector<Vector3> sites) {
    _edges.resize(state.edges.size());

    for (size_t k = 0; k < state.edges.size(); ++k) {
        _edges[k].reserve(state.edges[k].size());

        for (const EdgeState& edge : state.edges[k]) {
            const Vector3& neighbor = sites[edge.neighbor_index];
            _edges[k].push_back(EdgeGradient{
                edge.neighbor_index,
                sweep_rate(sites[k], neighbor, edge.integrals),
                sweep_rate(neighbor, sites[k], edge.integrals)
            });
        }
    }
}

inline std::vector<double> CapacityJacobian::apply(const std::vector<Vector3>& directions) const {
    std::vector<double> rates(_edges.size(), 0.0);

    for (size_t k = 0; k < _edges.size(); ++k) {
        for (const EdgeGradient& edge : _edges[k]) {
            rates[k] += edge.own.dot(directions[k]) - edge.neighbor.dot(directions[edge.neighbor_index]);
        }
    }

    return rates;
}

inline std::vector<Vector3> CapacityJacobian::transpose_apply(const std::vector<double>& values) const {
    std::vector<Vector3> gradients(_edges.size(), Vector3::Zero());

    for (size_t k = 0; k < _edges.size(); ++k) {
        for (const EdgeGradient& edge : _edges[k]) {
            gradients[k] += (values[k] - values[edge.neighbor_index]) * edge.own;
        }
    }

    return gradients;
}

// A displacement of the site moves the shared bisector outward at a rate
// proportional to the boundary point's projection on it. The site is
// subtracted so the rate is expressed about the site itself, which changes
// nothing for the tangential displacements the optimizer takes.
inline Vector3 CapacityJacobian::sweep_rate(
    const Vector3& site,
    const Vector3& neighbor,
    const fields::spherical::RegionIntegrals& integrals
) {
    double separation = (neighbor - site).norm();

    if (separation < GEOMETRIC_EPSILON) {
        return Vector3::Zero();
    }

    return (integrals.first_moment - site * integrals.mass) / separation;
}

} // namespace globe::voronoi::spherical

#endif //GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_CAPACITY_JACOBIAN_HPP_
