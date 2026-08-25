#ifndef GEOMETRY_ART_GEOMETRY_PLANAR_PERIODIC_DELAUNAY_HPP_
#define GEOMETRY_ART_GEOMETRY_PLANAR_PERIODIC_DELAUNAY_HPP_

#include "../../types.hpp"
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/assertions.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

namespace geometry_art::geometry::planar {

// One triangle of the fan around a site: the neighbour the spoke points at,
// and the circumcentre of the triangle lying between this spoke and the next.
// A site's Voronoi cell is exactly the loop of those circumcentres, and the
// bisector it shares with the neighbour runs between the circumcentres of the
// two spokes that bracket it.
struct Spoke {
    size_t neighbor_index;
    Vector2 neighbor_position;
    Vector2 circumcenter;
};

// A Delaunay triangulation of the flat torus: the plane modulo a rectangle
// of periods, in which every site has a finite fan of triangles and there is
// no boundary and no infinite face.
//
// The triangulation is of the periodically tiled sites, computed literally,
// with the cells read off the copies that lie in the domain. Interior sites
// need no copies at all -- a copy can only be a Delaunay neighbour of the
// domain's own sites if it lies within the largest circumradius of the
// domain, which is a few cell widths rather than a period -- so images are
// laid in a band around the seam whose width is measured rather than assumed.
//
// Positions come back in the querying site's own chart, so a neighbour a
// period away reads as the copy that actually borders it and no caller
// handles an offset.
class PeriodicDelaunay {
 public:
    PeriodicDelaunay(double width, double height, std::vector<Vector2> sites);

    PeriodicDelaunay(const PeriodicDelaunay&) = delete;
    PeriodicDelaunay& operator=(const PeriodicDelaunay&) = delete;
    PeriodicDelaunay(PeriodicDelaunay&&) = default;
    PeriodicDelaunay& operator=(PeriodicDelaunay&&) = default;

    [[nodiscard]] size_t size() const { return _handles.size(); }
    [[nodiscard]] double width() const { return _width; }
    [[nodiscard]] double height() const { return _height; }

    [[nodiscard]] Vector2 site(size_t index) const;
    [[nodiscard]] std::vector<Spoke> fan(size_t index) const;
    [[nodiscard]] Vector2 canonical(const Vector2& point) const;

 private:
    // A periodic image carries the index of the site it copies, so a wrapped
    // neighbour reports the same index its canonical placement would, and a
    // flag for the one copy that lies in the domain and owns the cell.
    struct VertexInfo {
        size_t site_index;
        bool central;
    };

    using Kernel = ::CGAL::Exact_predicates_inexact_constructions_kernel;
    using VertexBase = ::CGAL::Triangulation_vertex_base_with_info_2<VertexInfo, Kernel>;
    using DataStructure = ::CGAL::Triangulation_data_structure_2<VertexBase>;
    using Triangulation = ::CGAL::Delaunay_triangulation_2<Kernel, DataStructure>;
    using VertexHandle = Triangulation::Vertex_handle;
    using FaceHandle = Triangulation::Face_handle;

    // How many cell widths of images to lay around the domain before the
    // covering is measured. Wide enough that one pass almost always suffices,
    // and cheap enough that overshooting costs little.
    static constexpr double BAND_CELL_LAYERS = 4.0;

    double _width;
    double _height;
    std::unique_ptr<Triangulation> _triangulation;
    std::vector<VertexHandle> _handles;

    void build(const std::vector<Vector2>& sites, double band);
    [[nodiscard]] double initial_band(size_t count) const;
    [[nodiscard]] double distance_to_domain(const Vector2& point) const;
    [[nodiscard]] double maximum_circumradius() const;
    [[nodiscard]] double circumradius(FaceHandle face) const;
    [[nodiscard]] Vector2 position(VertexHandle vertex) const;
    [[nodiscard]] Vector2 dual_vertex(FaceHandle face) const;
};

// The band is laid at a guess, measured against the circumradii it produced,
// and widened until it covers them -- terminating at the full nine-tile
// replication, which is the most a covering can be and which suffices exactly
// when no circumradius reaches half a period.
inline PeriodicDelaunay::PeriodicDelaunay(double width, double height, std::vector<Vector2> sites) :
    _width(width),
    _height(height),
    _triangulation(std::make_unique<Triangulation>()) {
    CGAL_precondition(width > 0.0 && height > 0.0);

    for (Vector2& site : sites) {
        site = canonical(site);
    }

    double full_replication = std::max(_width, _height);
    double band = initial_band(sites.size());

    while (true) {
        build(sites, band);
        double radius = maximum_circumradius();

        if (radius < band || band >= full_replication) {
            CGAL_assertion(radius < 0.5 * std::min(_width, _height));
            static_cast<void>(radius);
            return;
        }

        band = std::min(2.0 * band, full_replication);
    }
}

inline Vector2 PeriodicDelaunay::site(size_t index) const {
    return position(_handles[index]);
}

inline std::vector<Spoke> PeriodicDelaunay::fan(size_t index) const {
    std::vector<Spoke> spokes;
    VertexHandle vertex = _handles[index];
    auto circulator = _triangulation->incident_faces(vertex);
    auto start = circulator;

    do {
        FaceHandle face = circulator;
        CGAL_precondition(!_triangulation->is_infinite(face));

        // The face's vertices run counter-clockwise, so the neighbour one
        // step round from the site opens the sector this face spans.
        VertexHandle neighbor = face->vertex(Triangulation::ccw(face->index(vertex)));

        spokes.push_back(Spoke{
            neighbor->info().site_index,
            position(neighbor),
            dual_vertex(face)
        });

        ++circulator;
    } while (circulator != start);

    return spokes;
}

inline Vector2 PeriodicDelaunay::canonical(const Vector2& point) const {
    double x = std::fmod(point.x(), _width);
    double y = std::fmod(point.y(), _height);

    if (x < 0.0) { x += _width; }
    if (y < 0.0) { y += _height; }

    return Vector2(x, y);
}

inline double PeriodicDelaunay::initial_band(size_t count) const {
    if (count == 0) {
        return std::max(_width, _height);
    }

    return BAND_CELL_LAYERS * std::sqrt(_width * _height / static_cast<double>(count));
}

inline void PeriodicDelaunay::build(const std::vector<Vector2>& sites, double band) {
    _triangulation = std::make_unique<Triangulation>();
    _handles.assign(sites.size(), VertexHandle());

    std::vector<std::pair<Triangulation::Point, VertexInfo>> points;
    points.reserve(sites.size());

    for (size_t index = 0; index < sites.size(); ++index) {
        for (int row = -1; row <= 1; ++row) {
            for (int column = -1; column <= 1; ++column) {
                bool central = row == 0 && column == 0;
                Vector2 copy = sites[index] + Vector2(column * _width, row * _height);

                if (!central && distance_to_domain(copy) > band) {
                    continue;
                }

                points.emplace_back(
                    Triangulation::Point(copy.x(), copy.y()),
                    VertexInfo{index, central}
                );
            }
        }
    }

    // The range insertion spatially sorts the points, so each location walk
    // starts beside the last instead of crossing the triangulation.
    _triangulation->insert(points.begin(), points.end());

    for (auto vertex = _triangulation->finite_vertices_begin();
        vertex != _triangulation->finite_vertices_end(); ++vertex) {
        if (vertex->info().central) {
            _handles[vertex->info().site_index] = vertex;
        }
    }
}

inline double PeriodicDelaunay::distance_to_domain(const Vector2& point) const {
    double across = std::max({0.0, -point.x(), point.x() - _width});
    double along = std::max({0.0, -point.y(), point.y() - _height});

    return std::sqrt(across * across + along * along);
}

// A site whose own fan touches the edge of the covered region has neighbours
// the band left out, which reads as an infinite incident face.
inline double PeriodicDelaunay::maximum_circumradius() const {
    double maximum = 0.0;

    for (const VertexHandle& vertex : _handles) {
        if (vertex == VertexHandle()) {
            continue;
        }

        auto circulator = _triangulation->incident_faces(vertex);
        auto start = circulator;

        do {
            if (_triangulation->is_infinite(circulator)) {
                return std::numeric_limits<double>::infinity();
            }

            maximum = std::max(maximum, circumradius(circulator));
            ++circulator;
        } while (circulator != start);
    }

    return maximum;
}

inline double PeriodicDelaunay::circumradius(FaceHandle face) const {
    auto center = _triangulation->dual(face);
    auto corner = _triangulation->point(face->vertex(0));

    return std::sqrt(::CGAL::to_double(::CGAL::squared_distance(center, corner)));
}

inline Vector2 PeriodicDelaunay::position(VertexHandle vertex) const {
    auto point = _triangulation->point(vertex);
    return Vector2(point.x(), point.y());
}

inline Vector2 PeriodicDelaunay::dual_vertex(FaceHandle face) const {
    auto center = _triangulation->dual(face);
    return Vector2(center.x(), center.y());
}

} // namespace geometry_art::geometry::planar

#endif //GEOMETRY_ART_GEOMETRY_PLANAR_PERIODIC_DELAUNAY_HPP_
