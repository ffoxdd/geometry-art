#ifndef GLOBEART_SRC_GLOBE_VORONOI_FLAT_CORE_TORUS_HPP_
#define GLOBEART_SRC_GLOBE_VORONOI_FLAT_CORE_TORUS_HPP_

#include "../../../types.hpp"
#include "../../../geometry/planar/polygon.hpp"
#include "../../../geometry/planar/segment.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/assertions.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

namespace globe::voronoi::flat {

using geometry::planar::Polygon;
using geometry::planar::Segment;

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

// A Voronoi diagram on the flat torus: the plane modulo a rectangle of
// periods. The diagram is the planar diagram of the periodically tiled
// sites, computed literally -- each site is inserted into the nine
// surrounding tiles and the cells are read off the central copy. That is
// exact as long as every Delaunay circumradius stays below half the
// shortest period, which the queries assert.
class Torus {
 public:
    Torus(double width, double height, std::vector<Vector2> sites);

    Torus(const Torus&) = delete;
    Torus& operator=(const Torus&) = delete;
    Torus(Torus&&) = default;
    Torus& operator=(Torus&&) = default;

    [[nodiscard]] size_t size() const { return _handles.size(); }
    [[nodiscard]] double width() const { return _width; }
    [[nodiscard]] double height() const { return _height; }
    [[nodiscard]] double area() const { return _width * _height; }

    [[nodiscard]] Vector2 site(size_t index) const;
    [[nodiscard]] Vector3 site_vector(size_t index) const;
    [[nodiscard]] std::vector<CellEdgeInfo> cell_edges(size_t index) const;
    [[nodiscard]] Polygon cell(size_t index) const;

    [[nodiscard]] Vector2 canonical(const Vector2& point) const;
    [[nodiscard]] std::unique_ptr<Torus> rebuilt(const std::vector<Vector3>& points) const;

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
    using Edge = std::pair<FaceHandle, int>;

    // How many cell widths of periodic images to lay around the domain before
    // the covering is measured. Wide enough that one pass almost always
    // suffices, and cheap enough that overshooting costs little.
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
    [[nodiscard]] Vector2 dual_vertex(FaceHandle face) const;
    [[nodiscard]] std::vector<Edge> incident_edges(VertexHandle vertex) const;
    void assert_covering_is_sufficient(FaceHandle face) const;
};

// Only sites near a seam need periodic images: a copy can be a Delaunay
// neighbour of the domain's own sites only if it lies within the largest
// circumradius of the domain, and that is a few cell widths, not a period.
// The band is therefore laid at a guess, measured against the circumradii it
// produced, and widened until it covers them -- terminating at the full
// nine-tile replication, which always suffices.
inline Torus::Torus(double width, double height, std::vector<Vector2> sites) :
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

        if (band >= full_replication || maximum_circumradius() < band) {
            return;
        }

        band = std::min(2.0 * band, full_replication);
    }
}

inline double Torus::initial_band(size_t count) const {
    if (count == 0) {
        return std::max(_width, _height);
    }

    return BAND_CELL_LAYERS * std::sqrt(area() / static_cast<double>(count));
}

inline void Torus::build(const std::vector<Vector2>& sites, double band) {
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

inline double Torus::distance_to_domain(const Vector2& point) const {
    double across = std::max({0.0, -point.x(), point.x() - _width});
    double along = std::max({0.0, -point.y(), point.y() - _height});

    return std::sqrt(across * across + along * along);
}

// A site whose own cell touches the edge of the covered region has neighbours
// the band left out, which reads as an infinite incident face.
inline double Torus::maximum_circumradius() const {
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

inline double Torus::circumradius(FaceHandle face) const {
    auto center = _triangulation->dual(face);
    auto corner = _triangulation->point(face->vertex(0));

    return std::sqrt(::CGAL::to_double(::CGAL::squared_distance(center, corner)));
}

inline Vector2 Torus::canonical(const Vector2& point) const {
    double x = std::fmod(point.x(), _width);
    double y = std::fmod(point.y(), _height);

    if (x < 0.0) { x += _width; }
    if (y < 0.0) { y += _height; }

    return Vector2(x, y);
}

inline Vector2 Torus::site(size_t index) const {
    auto point = _triangulation->point(_handles[index]);
    return Vector2(point.x(), point.y());
}

inline Vector3 Torus::site_vector(size_t index) const {
    Vector2 point = site(index);
    return Vector3(point.x(), point.y(), 0.0);
}

inline std::unique_ptr<Torus> Torus::rebuilt(const std::vector<Vector3>& points) const {
    std::vector<Vector2> sites;
    sites.reserve(points.size());

    for (const Vector3& point : points) {
        sites.emplace_back(point.x(), point.y());
    }

    return std::make_unique<Torus>(_width, _height, std::move(sites));
}

inline std::vector<Torus::Edge> Torus::incident_edges(VertexHandle vertex) const {
    std::vector<Edge> edges;
    auto circulator = _triangulation->incident_edges(vertex);
    auto start = circulator;

    do {
        edges.emplace_back(circulator->first, circulator->second);
        ++circulator;
    } while (circulator != start);

    return edges;
}

inline std::vector<CellEdgeInfo> Torus::cell_edges(size_t index) const {
    std::vector<CellEdgeInfo> result;
    VertexHandle vertex = _handles[index];

    for (const Edge& edge : incident_edges(vertex)) {
        FaceHandle face = edge.first;
        int position = edge.second;
        CGAL_precondition(!_triangulation->is_infinite(face));

        VertexHandle first = face->vertex(Triangulation::cw(position));
        VertexHandle second = face->vertex(Triangulation::ccw(position));
        VertexHandle neighbor = first == vertex ? second : first;

        FaceHandle target_face = face->neighbor(position);
        CGAL_precondition(!_triangulation->is_infinite(target_face));
        assert_covering_is_sufficient(face);

        auto neighbor_point = _triangulation->point(neighbor);
        auto source_opposite_point = _triangulation->point(face->vertex(position));
        auto target_opposite_point = _triangulation->point(target_face->vertex(target_face->index(face)));

        result.push_back(CellEdgeInfo{
            neighbor->info().site_index,
            face->vertex(position)->info().site_index,
            target_face->vertex(target_face->index(face))->info().site_index,
            Vector2(neighbor_point.x(), neighbor_point.y()),
            Vector2(source_opposite_point.x(), source_opposite_point.y()),
            Vector2(target_opposite_point.x(), target_opposite_point.y()),
            Segment(dual_vertex(face), dual_vertex(target_face))
        });
    }

    return result;
}

inline Polygon Torus::cell(size_t index) const {
    std::vector<Vector2> vertices;
    VertexHandle vertex = _handles[index];
    auto circulator = _triangulation->incident_faces(vertex);
    auto start = circulator;

    do {
        FaceHandle face = circulator;
        CGAL_precondition(!_triangulation->is_infinite(face));
        vertices.push_back(dual_vertex(face));
        ++circulator;
    } while (circulator != start);

    Polygon polygon(vertices);

    if (polygon.area() < 0.0) {
        std::reverse(vertices.begin(), vertices.end());
        return Polygon(std::move(vertices));
    }

    return polygon;
}

inline Vector2 Torus::dual_vertex(FaceHandle face) const {
    auto center = _triangulation->dual(face);
    return Vector2(center.x(), center.y());
}

inline void Torus::assert_covering_is_sufficient(FaceHandle face) const {
    CGAL_assertion(circumradius(face) < 0.5 * std::min(_width, _height));
    static_cast<void>(face);
}

} // namespace globe::voronoi::flat

#endif //GLOBEART_SRC_GLOBE_VORONOI_FLAT_CORE_TORUS_HPP_
