#ifndef GEOMETRY_ART_VORONOI_HESSIAN_BLOCKS_HPP_
#define GEOMETRY_ART_VORONOI_HESSIAN_BLOCKS_HPP_

#include "block_preconditioner.hpp"
#include "../types.hpp"
#include <CGAL/assertions.h>
#include <cstddef>
#include <utility>
#include <vector>

namespace geometry_art::voronoi {

struct NeighborBlock {
    size_t neighbor_index;
    Matrix3 value;
};

// A sparse symmetric operator over the sites as ambient vectors. Only
// Delaunay neighbours couple, so the matrix is stored as one dense block per
// cell plus one per shared bisector.
//
// through_manifold carries the blocks into the parameterisation the
// optimizer descends, restricted to each site's tangent space.
struct HessianBlocks {
    std::vector<Matrix3> diagonal;
    std::vector<std::vector<NeighborBlock>> neighbors;

    [[nodiscard]] std::vector<Vector3> multiply(const std::vector<Vector3>& directions) const;
    [[nodiscard]] HessianBlocks plus(const HessianBlocks& other) const;

    template<typename ManifoldType>
    [[nodiscard]] HessianBlocks through_manifold(
        const std::vector<Vector3>& points,
        const std::vector<Vector3>& site_gradients
    ) const;
};

inline std::vector<Vector3> HessianBlocks::multiply(const std::vector<Vector3>& directions) const {
    std::vector<Vector3> result(directions.size());

    for (size_t k = 0; k < directions.size(); ++k) {
        Vector3 sum = diagonal[k] * directions[k];

        for (const NeighborBlock& block : neighbors[k]) {
            sum += block.value * directions[block.neighbor_index];
        }

        result[k] = sum;
    }

    return result;
}

// Both operands must come from the same diagram, so the blocks line up
// position by position.
inline HessianBlocks HessianBlocks::plus(const HessianBlocks& other) const {
    CGAL_precondition(diagonal.size() == other.diagonal.size());

    HessianBlocks result = *this;

    for (size_t k = 0; k < diagonal.size(); ++k) {
        result.diagonal[k] += other.diagonal[k];
        CGAL_precondition(neighbors[k].size() == other.neighbors[k].size());

        for (size_t position = 0; position < neighbors[k].size(); ++position) {
            CGAL_precondition(
                neighbors[k][position].neighbor_index == other.neighbors[k][position].neighbor_index
            );
            result.neighbors[k][position].value += other.neighbors[k][position].value;
        }
    }

    return result;
}

template<typename ManifoldType>
HessianBlocks HessianBlocks::through_manifold(
    const std::vector<Vector3>& points,
    const std::vector<Vector3>& site_gradients
) const {
    std::vector<ManifoldType> manifolds;
    manifolds.reserve(points.size());

    for (const Vector3& point : points) {
        manifolds.emplace_back(point);
    }

    HessianBlocks result;
    result.diagonal.reserve(points.size());
    result.neighbors.resize(points.size());

    for (size_t k = 0; k < points.size(); ++k) {
        result.diagonal.push_back(manifolds[k].tangential_hessian(site_gradients[k], diagonal[k]));
        result.neighbors[k].reserve(neighbors[k].size());

        for (const NeighborBlock& block : neighbors[k]) {
            result.neighbors[k].push_back(NeighborBlock{
                block.neighbor_index,
                manifolds[k].mixed_hessian(block.value, manifolds[block.neighbor_index])
            });
        }
    }

    return result;
}

// The blocks paired with the inverse of their block diagonal, which is what
// the descent's conjugate gradient asks a curvature operator for. The pairing
// is explicit because the inverse is built once the blocks are final, not on
// every product.
struct PreconditionedBlocks {
    HessianBlocks blocks;
    BlockPreconditioner preconditioner;

    [[nodiscard]] static PreconditionedBlocks of(HessianBlocks blocks);

    [[nodiscard]] std::vector<Vector3> multiply(const std::vector<Vector3>& directions) const {
        return blocks.multiply(directions);
    }

    [[nodiscard]] std::vector<Vector3> precondition(const std::vector<Vector3>& residuals) const {
        return preconditioner.apply(residuals);
    }
};

inline PreconditionedBlocks PreconditionedBlocks::of(HessianBlocks blocks) {
    BlockPreconditioner preconditioner(blocks.diagonal);
    return PreconditionedBlocks{std::move(blocks), std::move(preconditioner)};
}

} // namespace geometry_art::voronoi

#endif //GEOMETRY_ART_VORONOI_HESSIAN_BLOCKS_HPP_
