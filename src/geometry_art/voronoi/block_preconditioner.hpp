#ifndef GEOMETRY_ART_VORONOI_BLOCK_PRECONDITIONER_HPP_
#define GEOMETRY_ART_VORONOI_BLOCK_PRECONDITIONER_HPP_

#include "../types.hpp"
#include <Eigen/Eigenvalues>
#include <cmath>
#include <cstddef>
#include <vector>

namespace geometry_art::voronoi {

// The block-diagonal part of a curvature operator, inverted, which is what
// the conjugate gradient needs to stop resolving the spread in scale between
// sites instead of the shape of the problem.
//
// A block is inverted through its eigenvalues rather than directly, because
// two things have to be told apart. A site's step is confined to the
// parameterisation the geometry descends -- a tangent plane on the sphere,
// the plane itself on a flat domain -- so every block is rank deficient in
// the ambient third direction, and that direction is annihilated, since a
// true inverse would amplify it without bound. A direction of genuinely
// negative curvature is not absent, only inverted, so it is scaled by its
// magnitude and left in; dropping it would freeze the site until the
// curvature turned, and the descent has its own answer to negative
// curvature.
class BlockPreconditioner {
 public:
    BlockPreconditioner() = default;
    explicit BlockPreconditioner(const std::vector<Matrix3>& blocks);

    [[nodiscard]] std::vector<Vector3> apply(const std::vector<Vector3>& residuals) const;

 private:
    static constexpr double RELATIVE_EIGENVALUE_FLOOR = 1e-9;

    std::vector<Matrix3> _inverses;

    [[nodiscard]] static Matrix3 pseudo_inverse(const Matrix3& block);
};

inline BlockPreconditioner::BlockPreconditioner(const std::vector<Matrix3>& blocks) {
    _inverses.reserve(blocks.size());

    for (const Matrix3& block : blocks) {
        _inverses.push_back(pseudo_inverse(block));
    }
}

inline std::vector<Vector3> BlockPreconditioner::apply(const std::vector<Vector3>& residuals) const {
    if (_inverses.empty()) {
        return residuals;
    }

    std::vector<Vector3> result(residuals.size());

    for (size_t k = 0; k < residuals.size(); ++k) {
        result[k] = _inverses[k] * residuals[k];
    }

    return result;
}

inline Matrix3 BlockPreconditioner::pseudo_inverse(const Matrix3& block) {
    Matrix3 symmetric = 0.5 * (block + block.transpose());
    Eigen::SelfAdjointEigenSolver<Matrix3> solver(symmetric);

    if (solver.info() != Eigen::Success) {
        return Matrix3::Identity();
    }

    Vector3 eigenvalues = solver.eigenvalues();
    double largest = eigenvalues.cwiseAbs().maxCoeff();

    if (largest <= 0.0) {
        return Matrix3::Identity();
    }

    double floor = RELATIVE_EIGENVALUE_FLOOR * largest;
    Vector3 inverted = Vector3::Zero();

    for (int axis = 0; axis < 3; ++axis) {
        double magnitude = std::abs(eigenvalues[axis]);
        inverted[axis] = magnitude > floor ? 1.0 / magnitude : 0.0;
    }

    return solver.eigenvectors() * inverted.asDiagonal() * solver.eigenvectors().transpose();
}

} // namespace geometry_art::voronoi

#endif //GEOMETRY_ART_VORONOI_BLOCK_PRECONDITIONER_HPP_
