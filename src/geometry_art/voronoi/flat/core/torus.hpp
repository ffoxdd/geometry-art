#ifndef GEOMETRY_ART_VORONOI_FLAT_CORE_TORUS_HPP_
#define GEOMETRY_ART_VORONOI_FLAT_CORE_TORUS_HPP_

#include "../../../types.hpp"
#include "../../../geometry/planar/periodic_delaunay.hpp"
#include "../../../geometry/planar/polygon.hpp"
#include "../../../geometry/planar/segment.hpp"
#include <algorithm>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

namespace geometry_art::voronoi::flat {

using geometry::planar::PeriodicDelaunay;
using geometry::planar::Polygon;
using geometry::planar::Segment;
using geometry::planar::Spoke;

// The segment's endpoints are Voronoi vertices, each equidistant from the
// two sites of the edge and one more; the opposite indices name that third
// site. Everything geometric is expressed in the cell's own chart -- the
// plane around the cell's canonical site -- so the neighbor's position is
// carried here rather than read from its own canonical placement, which may
// sit a period away.
struct CellEdgeInfo {
    size_t neighbor_index;
    size_t source_opposite_index;
    size_t target_opposite_index;
    Vector2 neighbor_position;
    Vector2 source_opposite_position;
    Vector2 target_opposite_position;
    Segment boundary;
};

// A Voronoi diagram on the flat torus, read off the dual of a periodic
// Delaunay triangulation: a cell is the loop of circumcentres around its
// site, and the bisector it shares with a neighbour spans the two
// circumcentres that bracket that neighbour.
class Torus {
 public:
    Torus(double width, double height, std::vector<Vector2> sites);

    Torus(const Torus&) = delete;
    Torus& operator=(const Torus&) = delete;
    Torus(Torus&&) = default;
    Torus& operator=(Torus&&) = default;

    [[nodiscard]] size_t size() const { return _triangulation.size(); }
    [[nodiscard]] double width() const { return _triangulation.width(); }
    [[nodiscard]] double height() const { return _triangulation.height(); }
    [[nodiscard]] double area() const { return width() * height(); }

    [[nodiscard]] Vector2 site(size_t index) const { return _triangulation.site(index); }
    [[nodiscard]] Vector3 site_vector(size_t index) const;
    [[nodiscard]] std::vector<CellEdgeInfo> cell_edges(size_t index) const;
    [[nodiscard]] Polygon cell(size_t index) const;

    [[nodiscard]] Vector2 canonical(const Vector2& point) const {
        return _triangulation.canonical(point);
    }

    [[nodiscard]] std::unique_ptr<Torus> rebuilt(const std::vector<Vector3>& points) const;

 private:
    PeriodicDelaunay _triangulation;
};

inline Torus::Torus(double width, double height, std::vector<Vector2> sites) :
    _triangulation(width, height, std::move(sites)) {
}

inline Vector3 Torus::site_vector(size_t index) const {
    Vector2 point = site(index);
    return Vector3(point.x(), point.y(), 0.0);
}

inline std::vector<CellEdgeInfo> Torus::cell_edges(size_t index) const {
    std::vector<Spoke> spokes = _triangulation.fan(index);
    std::vector<CellEdgeInfo> result;
    result.reserve(spokes.size());

    size_t count = spokes.size();

    for (size_t position = 0; position < count; ++position) {
        const Spoke& previous = spokes[(position + count - 1) % count];
        const Spoke& spoke = spokes[position];
        const Spoke& next = spokes[(position + 1) % count];

        result.push_back(CellEdgeInfo{
            spoke.neighbor_index,
            previous.neighbor_index,
            next.neighbor_index,
            spoke.neighbor_position,
            previous.neighbor_position,
            next.neighbor_position,
            Segment(previous.circumcenter, spoke.circumcenter)
        });
    }

    return result;
}

inline Polygon Torus::cell(size_t index) const {
    std::vector<Vector2> vertices;

    for (const Spoke& spoke : _triangulation.fan(index)) {
        vertices.push_back(spoke.circumcenter);
    }

    Polygon polygon(vertices);

    if (polygon.area() < 0.0) {
        std::reverse(vertices.begin(), vertices.end());
        return Polygon(std::move(vertices));
    }

    return polygon;
}

inline std::unique_ptr<Torus> Torus::rebuilt(const std::vector<Vector3>& points) const {
    std::vector<Vector2> sites;
    sites.reserve(points.size());

    for (const Vector3& point : points) {
        sites.emplace_back(point.x(), point.y());
    }

    return std::make_unique<Torus>(width(), height(), std::move(sites));
}

} // namespace geometry_art::voronoi::flat

#endif //GEOMETRY_ART_VORONOI_FLAT_CORE_TORUS_HPP_
