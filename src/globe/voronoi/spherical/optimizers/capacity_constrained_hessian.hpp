#ifndef GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_CAPACITY_CONSTRAINED_HESSIAN_HPP_
#define GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_CAPACITY_CONSTRAINED_HESSIAN_HPP_

#include "../../block_preconditioner.hpp"
#include "../../capacity_jacobian.hpp"
#include "../../hessian_blocks.hpp"
#include "../../../types.hpp"
#include <cstddef>
#include <utility>
#include <vector>

namespace globe::voronoi::spherical {

// Curvature of the augmented Lagrangian, in the Gauss-Newton form: the exact
// CVT curvature plus the penalty's outer product of constraint gradients.
//
// The standard form would add the constraints' own second derivatives,
// weighted by the multipliers and the violations. Those need the velocities
// of the Voronoi vertices, which the moment machinery does not yet supply;
// near a solution the violations vanish while the penalty term is the one
// that grows, so the term kept is the one that governs the conditioning.
class CapacityConstrainedHessian {
 public:
    CapacityConstrainedHessian(
        HessianBlocks energy_curvature,
        CapacityJacobian jacobian,
        std::vector<Vector3> sites,
        double penalty
    );

    [[nodiscard]] std::vector<Vector3> multiply(const std::vector<Vector3>& directions) const;

    [[nodiscard]] std::vector<Vector3> precondition(const std::vector<Vector3>& residuals) const {
        return _preconditioner.apply(residuals);
    }

 private:
    HessianBlocks _energy_curvature;
    CapacityJacobian _jacobian;
    std::vector<Vector3> _sites;
    double _penalty;
    BlockPreconditioner _preconditioner;

    [[nodiscard]] std::vector<Vector3> tangential(const std::vector<Vector3>& value) const;
    [[nodiscard]] std::vector<Matrix3> diagonal_blocks() const;
};

inline CapacityConstrainedHessian::CapacityConstrainedHessian(
    HessianBlocks energy_curvature,
    CapacityJacobian jacobian,
    std::vector<Vector3> sites,
    double penalty
) :
    _energy_curvature(std::move(energy_curvature)),
    _jacobian(std::move(jacobian)),
    _sites(std::move(sites)),
    _penalty(penalty) {
    _preconditioner = BlockPreconditioner(diagonal_blocks());
}

inline std::vector<Vector3> CapacityConstrainedHessian::multiply(const std::vector<Vector3>& directions) const {
    std::vector<Vector3> tangent = tangential(directions);
    std::vector<Vector3> result = _energy_curvature.multiply(tangent);

    if (_penalty <= 0.0) {
        return result;
    }

    std::vector<Vector3> penalty_term = tangential(_jacobian.transpose_apply(_jacobian.apply(tangent)));

    for (size_t k = 0; k < result.size(); ++k) {
        result[k] += _penalty * penalty_term[k];
    }

    return result;
}

// The sites are directions, so anything along one of them is unobservable.
inline std::vector<Vector3> CapacityConstrainedHessian::tangential(const std::vector<Vector3>& value) const {
    std::vector<Vector3> result(value.size());

    for (size_t k = 0; k < value.size(); ++k) {
        result[k] = value[k] - value[k].dot(_sites[k]) * _sites[k];
    }

    return result;
}

// The penalty pulls the block diagonal apart as it grows, which is exactly
// what the conjugate gradient struggles with, so the diagonal it is
// preconditioned by carries the penalty's own share, projected onto the
// tangent planes the descent moves in.
inline std::vector<Matrix3> CapacityConstrainedHessian::diagonal_blocks() const {
    std::vector<Matrix3> result = _energy_curvature.diagonal;

    if (_penalty > 0.0) {
        std::vector<Matrix3> constraint = _jacobian.normal_equations_diagonal();

        for (size_t k = 0; k < result.size(); ++k) {
            result[k] += _penalty * constraint[k];
        }
    }

    for (size_t k = 0; k < result.size(); ++k) {
        Matrix3 projection = Matrix3::Identity() - _sites[k] * _sites[k].transpose();
        result[k] = projection * result[k] * projection;
    }

    return result;
}

} // namespace globe::voronoi::spherical

#endif //GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_CAPACITY_CONSTRAINED_HESSIAN_HPP_
