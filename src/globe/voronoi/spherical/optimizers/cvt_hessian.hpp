#ifndef GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_CVT_HESSIAN_HPP_
#define GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_CVT_HESSIAN_HPP_

#include "../../../types.hpp"
#include "../../../fields/spherical/field.hpp"
#include "../../../fields/spherical/polynomial_field.hpp"
#include "../../../math/normalization.hpp"
#include "../../../std_ext/parallel_for.hpp"
#include "../core/sphere.hpp"
#include <cstddef>
#include <utility>
#include <vector>

namespace globe::voronoi::spherical {

struct NeighborBlock {
    size_t neighbor_index;
    Matrix3 value;
};

// Second derivatives of the CVT energy with respect to the sites as ambient
// vectors. Only Delaunay neighbours couple, so the matrix is stored as one
// dense block per cell plus one per shared bisector.
//
// through_normalization carries the blocks into the parameterisation the
// optimizer descends, restricted to each site's tangent plane.
struct CvtHessianBlocks {
    std::vector<Matrix3> diagonal;
    std::vector<std::vector<NeighborBlock>> neighbors;

    [[nodiscard]] std::vector<Vector3> multiply(const std::vector<Vector3>& directions) const;
    [[nodiscard]] CvtHessianBlocks through_normalization(
        const std::vector<Vector3>& points,
        const std::vector<Vector3>& site_gradients
    ) const;
};

template<fields::spherical::Field FieldType = fields::spherical::PolynomialField>
class CvtHessian {
 public:
    explicit CvtHessian(FieldType field);

    [[nodiscard]] CvtHessianBlocks assemble(const Sphere& sphere) const;

 private:
    FieldType _field;

    [[nodiscard]] Matrix3 bisector_block(const Vector3& site, const Vector3& neighbor, const Arc& arc) const;
};

inline std::vector<Vector3> CvtHessianBlocks::multiply(const std::vector<Vector3>& directions) const {
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

inline CvtHessianBlocks CvtHessianBlocks::through_normalization(
    const std::vector<Vector3>& points,
    const std::vector<Vector3>& site_gradients
) const {
    std::vector<Normalization> normalizations;
    normalizations.reserve(points.size());

    for (const Vector3& point : points) {
        normalizations.emplace_back(point);
    }

    CvtHessianBlocks result;
    result.diagonal.reserve(points.size());
    result.neighbors.resize(points.size());

    for (size_t k = 0; k < points.size(); ++k) {
        result.diagonal.push_back(normalizations[k].tangential_hessian(site_gradients[k], diagonal[k]));
        result.neighbors[k].reserve(neighbors[k].size());

        for (const NeighborBlock& block : neighbors[k]) {
            result.neighbors[k].push_back(NeighborBlock{
                block.neighbor_index,
                normalizations[k].mixed_hessian(block.value, normalizations[block.neighbor_index])
            });
        }
    }

    return result;
}

template<fields::spherical::Field FieldType>
CvtHessian<FieldType>::CvtHessian(FieldType field) :
    _field(std::move(field)) {
}

template<fields::spherical::Field FieldType>
CvtHessianBlocks CvtHessian<FieldType>::assemble(const Sphere& sphere) const {
    size_t count = sphere.size();
    std::vector<Vector3> sites(count);
    std::vector<std::vector<CellEdgeInfo>> cell_edges(count);

    for (size_t k = 0; k < count; ++k) {
        sites[k] = to_vector3(sphere.site(k));
        cell_edges[k] = sphere.cell_edges(k);
    }

    CvtHessianBlocks blocks;
    blocks.diagonal.assign(count, Matrix3::Zero());
    blocks.neighbors.resize(count);

    std_ext::parallel_for(count, [&](size_t k) {
        blocks.neighbors[k].reserve(cell_edges[k].size());

        for (const CellEdgeInfo& edge : cell_edges[k]) {
            Matrix3 block = 2.0 * bisector_block(sites[k], sites[edge.neighbor_index], edge.arc);
            blocks.diagonal[k] -= block;
            blocks.neighbors[k].push_back(NeighborBlock{edge.neighbor_index, block});
        }
    });

    return blocks;
}

// Moving a site sweeps the shared bisector at a rate proportional to the
// point's projection on the displacement, so the block is the density's
// second moment along the arc scaled by the site separation.
template<fields::spherical::Field FieldType>
Matrix3 CvtHessian<FieldType>::bisector_block(
    const Vector3& site,
    const Vector3& neighbor,
    const Arc& arc
) const {
    double separation = (neighbor - site).norm();

    if (separation < GEOMETRIC_EPSILON) {
        return Matrix3::Zero();
    }

    return _field.second_moment(arc) / separation;
}

} // namespace globe::voronoi::spherical

#endif //GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_CVT_HESSIAN_HPP_
