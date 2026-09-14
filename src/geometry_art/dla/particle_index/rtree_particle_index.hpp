#ifndef GEOMETRY_ART_DLA_PARTICLE_INDEX_RTREE_PARTICLE_INDEX_HPP_
#define GEOMETRY_ART_DLA_PARTICLE_INDEX_RTREE_PARTICLE_INDEX_HPP_

#include "particle_index.hpp"
#include "../../types.hpp"
#include <CGAL/assertions.h>
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <cstddef>
#include <utility>

namespace geometry_art::dla {

// An R-tree rather than the kd-tree used elsewhere in the library: CGAL's
// Kd_tree rebuilds on the first query after an insertion, which the
// one-particle-at-a-time growth here would trigger on every step.
class RTreeParticleIndex {
 public:
    void insert(const Vector3& center, std::size_t index);
    [[nodiscard]] NearestParticle nearest(const Vector3& point) const;

 private:
    using Point = boost::geometry::model::point<double, 3, boost::geometry::cs::cartesian>;
    using Entry = std::pair<Point, std::size_t>;
    using Tree = boost::geometry::index::rtree<Entry, boost::geometry::index::quadratic<16>>;

    Tree _tree;

    [[nodiscard]] static Point to_point(const Vector3& vector);
    [[nodiscard]] static Vector3 to_vector(const Point& point);
};

inline void RTreeParticleIndex::insert(const Vector3& center, std::size_t index) {
    _tree.insert(Entry(to_point(center), index));
}

inline NearestParticle RTreeParticleIndex::nearest(const Vector3& point) const {
    CGAL_precondition(!_tree.empty());

    Entry closest = *_tree.qbegin(boost::geometry::index::nearest(to_point(point), 1));

    return NearestParticle{closest.second, (point - to_vector(closest.first)).norm()};
}

inline RTreeParticleIndex::Point RTreeParticleIndex::to_point(const Vector3& vector) {
    return Point(vector.x(), vector.y(), vector.z());
}

inline Vector3 RTreeParticleIndex::to_vector(const Point& point) {
    return Vector3(
        boost::geometry::get<0>(point),
        boost::geometry::get<1>(point),
        boost::geometry::get<2>(point)
    );
}

} // namespace geometry_art::dla

#endif //GEOMETRY_ART_DLA_PARTICLE_INDEX_RTREE_PARTICLE_INDEX_HPP_
