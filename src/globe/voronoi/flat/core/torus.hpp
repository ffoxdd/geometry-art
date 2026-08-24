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
#include <memory>
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
    Torus(double width, double height);

    Torus(const Torus&) = delete;
    Torus& operator=(const Torus&) = delete;
    Torus(Torus&&) = default;
    Torus& operator=(Torus&&) = default;

    void insert(const Vector2& point);

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
    using Kernel = ::CGAL::Exact_predicates_inexact_constructions_kernel;
    using VertexBase = ::CGAL::Triangulation_vertex_base_with_info_2<size_t, Kernel>;
    using DataStructure = ::CGAL::Triangulation_data_structure_2<VertexBase>;
    using Triangulation = ::CGAL::Delaunay_triangulation_2<Kernel, DataStructure>;
    using VertexHandle = Triangulation::Vertex_handle;
    using FaceHandle = Triangulation::Face_handle;
    using Edge = std::pair<FaceHandle, int>;

    double _width;
    double _height;
    std::unique_ptr<Triangulation> _triangulation;
    std::vector<VertexHandle> _handles;

    [[nodiscard]] Vector2 dual_vertex(FaceHandle face) const;
    [[nodiscard]] std::vector<Edge> incident_edges(VertexHandle vertex) const;
    void assert_covering_is_sufficient(FaceHandle face) const;
};

inline Torus::Torus(double width, double height) :
    _width(width),
    _height(height),
    _triangulation(std::make_unique<Triangulation>()) {
    CGAL_precondition(width > 0.0 && height > 0.0);
}

inline void Torus::insert(const Vector2& point) {
    Vector2 base = canonical(point);
    size_t index = _handles.size();
    VertexHandle center = nullptr;

    for (int row = -1; row <= 1; ++row) {
        for (int column = -1; column <= 1; ++column) {
            Vector2 copy = base + Vector2(column * _width, row * _height);
            VertexHandle handle = _triangulation->insert(Triangulation::Point(copy.x(), copy.y()));
            handle->info() = index;

            if (row == 0 && column == 0) {
                center = handle;
            }
        }
    }

    _handles.push_back(center);
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
    auto torus = std::make_unique<Torus>(_width, _height);

    for (const Vector3& point : points) {
        torus->insert(Vector2(point.x(), point.y()));
    }

    return torus;
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

        result.push_back(CellEdgeInfo{
            neighbor->info(),
            face->vertex(position)->info(),
            target_face->vertex(target_face->index(face))->info(),
            Vector2(neighbor_point.x(), neighbor_point.y()),
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
    auto center = _triangulation->dual(face);
    auto corner = _triangulation->point(face->vertex(0));
    double radius = std::sqrt(::CGAL::to_double(::CGAL::squared_distance(center, corner)));
    CGAL_assertion(radius < 0.5 * std::min(_width, _height));
    (void)radius;
}

} // namespace globe::voronoi::flat

#endif //GLOBEART_SRC_GLOBE_VORONOI_FLAT_CORE_TORUS_HPP_
