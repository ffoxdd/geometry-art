#ifndef GLOBEART_SRC_GLOBE_VORONOI_CAPACITY_JACOBIAN_HPP_
#define GLOBEART_SRC_GLOBE_VORONOI_CAPACITY_JACOBIAN_HPP_

#include "state.hpp"
#include "../types.hpp"
#include "../fields/region_integrals.hpp"
#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

namespace globe::voronoi {

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
    explicit CapacityJacobian(const DiagramState& state);

    [[nodiscard]] std::vector<double> apply(const std::vector<Vector3>& directions) const;
    [[nodiscard]] std::vector<Vector3> transpose_apply(const std::vector<double>& values) const;
    [[nodiscard]] std::vector<Matrix3> normal_equations_diagonal() const;

 private:
    struct EdgeGradient {
        size_t neighbor_index;
        Vector3 own;
        Vector3 neighbor;
    };

    std::vector<std::vector<EdgeGradient>> _edges;
};

// A displacement of a site moves the shared bisector outward at a normal
// speed of (x - site) . displacement / separation, the same expression on
// the sphere and on any flat domain. The state stores exactly the moments
// about each site, so no site coordinate enters here.
inline CapacityJacobian::CapacityJacobian(const DiagramState& state) {
    _edges.resize(state.edges.size());

    for (size_t k = 0; k < state.edges.size(); ++k) {
        _edges[k].reserve(state.edges[k].size());

        for (const EdgeState& edge : state.edges[k]) {
            if (edge.separation < GEOMETRIC_EPSILON) {
                _edges[k].push_back(EdgeGradient{edge.neighbor_index, Vector3::Zero(), Vector3::Zero()});
                continue;
            }

            _edges[k].push_back(EdgeGradient{
                edge.neighbor_index,
                edge.moment_about_own / edge.separation,
                edge.moment_about_neighbor / edge.separation
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

// The diagonal blocks of J^T J, which is the curvature the penalty term
// contributes. A cell's mass responds to its own site through the sum of its
// edge rates and to a neighbour's through that edge's rate, so each block is
// the outer square of one such sum. The sums are accumulated per site pair
// first: on a wrapped domain a cell can border the same site more than once,
// and the two edges answer to one column of the Jacobian.
inline std::vector<Matrix3> CapacityJacobian::normal_equations_diagonal() const {
    std::vector<Matrix3> result(_edges.size(), Matrix3::Zero());
    std::vector<std::pair<size_t, Vector3>> sums;

    for (size_t k = 0; k < _edges.size(); ++k) {
        Vector3 own = Vector3::Zero();
        sums.clear();

        for (const EdgeGradient& edge : _edges[k]) {
            own += edge.own;
            auto entry = std::find_if(sums.begin(), sums.end(), [&](const auto& candidate) {
                return candidate.first == edge.neighbor_index;
            });

            if (entry == sums.end()) {
                sums.emplace_back(edge.neighbor_index, edge.neighbor);
                continue;
            }

            entry->second += edge.neighbor;
        }

        result[k] += own * own.transpose();

        for (const auto& [index, sum] : sums) {
            result[index] += sum * sum.transpose();
        }
    }

    return result;
}

} // namespace globe::voronoi

#endif //GLOBEART_SRC_GLOBE_VORONOI_CAPACITY_JACOBIAN_HPP_
