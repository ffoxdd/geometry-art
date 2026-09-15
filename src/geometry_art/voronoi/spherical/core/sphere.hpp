#ifndef GEOMETRY_ART_VORONOI_SPHERICAL_CORE_SPHERE_HPP_
#define GEOMETRY_ART_VORONOI_SPHERICAL_CORE_SPHERE_HPP_

#include "../../../cgal/types.hpp"
#include "../../../geometry/spherical/arc.hpp"
#include "../../../geometry/spherical/cap.hpp"
#include "../../../geometry/spherical/cap_polygon.hpp"
#include <CGAL/Exact_spherical_kernel_3.h>
#include <CGAL/Delaunay_triangulation_on_sphere_traits_2.h>
#include <CGAL/Delaunay_triangulation_on_sphere_2.h>
#include "../../../geometry/spherical/polygon/polygon.hpp"
#include "circulator_iterator.hpp"
#include <algorithm>
#include <cstddef>
#include <memory>
#include <ranges>
#include <set>
#include <unordered_map>

namespace geometry_art::voronoi::spherical {

using geometry::spherical::Arc;
using geometry::spherical::Cap;
using geometry::spherical::CapPolygon;
using geometry::spherical::Polygon;

// The arc's endpoints are Voronoi vertices, each equidistant from the two
// sites of the edge and one more; the opposite indices name that third
// site, which is what the endpoint's motion depends on.
struct CellEdgeInfo {
    size_t neighbor_index;
    size_t source_opposite_index;
    size_t target_opposite_index;
    Arc arc;
};

struct VoronoiVertex {
    VectorS2 position;
    std::vector<Arc> arcs;  // arcs radiating from this vertex, ordered
};

class Sphere {
 public:
    explicit Sphere();

    Sphere(const Sphere&) = delete;
    Sphere& operator=(const Sphere&) = delete;
    Sphere(Sphere&&) = default;
    Sphere& operator=(Sphere&&) = default;

    void insert(cgal::Point3 point);
    std::size_t size() const;
    std::size_t cell_count() const;

    Polygon cell(size_t index) const;
    CapPolygon offset_cell(size_t index, double inset) const;
    auto cells() const;
    auto arcs() const;
    std::vector<VoronoiVertex> vertices() const;
    std::vector<Arc> unique_arcs() const;

    cgal::Point3 site(size_t index) const;
    void update_site(size_t index, cgal::Point3 new_position);

    std::vector<CellEdgeInfo> cell_edges(size_t index) const;

 private:
    using Kernel = ::CGAL::Exact_predicates_inexact_constructions_kernel;
    using SphericalKernel = ::CGAL::Exact_spherical_kernel_3;
    using Triangulation = ::CGAL::Delaunay_triangulation_on_sphere_2<
        ::CGAL::Delaunay_triangulation_on_sphere_traits_2<Kernel, SphericalKernel>
    >;
    using VertexHandle = Triangulation::Vertex_handle;
    using FaceHandle = Triangulation::Face_handle;
    using EdgeCirculator = Triangulation::Edge_circulator;
    using Edge = Triangulation::Edge;
    using EdgeCirculatorIterator = CirculatorIterator<EdgeCirculator, Edge>;

    std::unique_ptr<Triangulation> _triangulation;
    std::vector<VertexHandle> _handles;
    std::unordered_map<VertexHandle, size_t> _handle_to_index;

    auto static edge_circulator_range(EdgeCirculator edge_circulator);
    auto incident_edges_range(VertexHandle vertex_handle) const;
    std::vector<Arc> cell_arcs(size_t index) const;
    size_t vertex_index(VertexHandle handle) const;

    [[nodiscard]] VectorS2 dual_vertex(FaceHandle face) const;
    [[nodiscard]] Arc dual_arc(const Edge& edge) const;
};

inline Sphere::Sphere() :
    _triangulation(std::make_unique<Triangulation>()) {
}

inline void Sphere::insert(cgal::Point3 point) {
    VertexHandle handle = _triangulation->insert(point);
    size_t index = _handles.size();
    _handles.push_back(handle);
    _handle_to_index[handle] = index;
}

inline auto Sphere::edge_circulator_range(EdgeCirculator edge_circulator) {
    auto begin = EdgeCirculatorIterator(edge_circulator);
    auto end = EdgeCirculatorIterator(edge_circulator);

    return std::ranges::subrange(begin, end);
}

inline auto Sphere::incident_edges_range(VertexHandle vertex_handle) const {
    return edge_circulator_range(_triangulation->incident_edges(vertex_handle));
}

inline std::size_t Sphere::size() const {
    return _triangulation->number_of_vertices();
}

// Sites have cells once the triangulation spans the sphere.
inline std::size_t Sphere::cell_count() const {
    return _triangulation->dimension() >= 2 ? size() : 0;
}

inline cgal::Point3 Sphere::site(size_t index) const {
    return _triangulation->point(_handles[index]);
}

inline void Sphere::update_site(size_t index, cgal::Point3 new_position) {
    _handle_to_index.erase(_handles[index]);
    _triangulation->remove(_handles[index]);
    _handles[index] = _triangulation->insert(new_position);
    _handle_to_index[_handles[index]] = index;
}

inline std::vector<Arc> Sphere::cell_arcs(size_t index) const {
    std::vector<Arc> arcs;

    for (const CellEdgeInfo& edge : cell_edges(index)) {
        arcs.push_back(edge.arc);
    }

    if (arcs.size() < 2) {
        return arcs;
    }

    auto min_it = std::min_element(arcs.begin(), arcs.end(),
        [](const Arc& a, const Arc& b) {
            const auto& sa = a.source();
            const auto& sb = b.source();
            if (sa.x() != sb.x()) return sa.x() < sb.x();
            if (sa.y() != sb.y()) return sa.y() < sb.y();
            return sa.z() < sb.z();
        });

    std::rotate(arcs.begin(), min_it, arcs.end());

    return arcs;
}

inline Polygon Sphere::cell(size_t index) const {
    return Polygon(cell_arcs(index));
}

// The cell shrunk by a geodesic distance from every bisector: the
// intersection of its neighbours' hemispheres, each inset by that distance,
// whose rims are small circles rather than the bisectors' great circles.
inline CapPolygon Sphere::offset_cell(size_t index, double inset) const {
    if (cell_count() == 0) {
        return CapPolygon(std::vector<CapPolygon::Edge>{});
    }

    CapPolygon region = CapPolygon::of(cell(index));
    VectorS2 own = to_vector_s2(site(index));

    for (const CellEdgeInfo& edge : cell_edges(index)) {
        VectorS2 neighbor = to_vector_s2(site(edge.neighbor_index));
        region.clip(Cap::hemisphere_inset_by((own - neighbor).normalized(), inset));
    }

    return region;
}

inline auto Sphere::cells() const {
    return std::views::iota(size_t(0), cell_count()) | std::views::transform(
        [this](size_t index) {
            return cell(index);
        }
    );
}

inline auto Sphere::arcs() const {
    return cells() | std::views::transform(
        [](const Polygon &cell) {
            return cell.arcs();
        }
    ) | std::views::join;
}

inline std::vector<VoronoiVertex> Sphere::vertices() const {
    if (_triangulation->dimension() < 2) {
        return {};
    }

    std::vector<VoronoiVertex> result;

    // Voronoi vertices are duals of Delaunay faces
    for (auto fit = _triangulation->solid_faces_begin();
         fit != _triangulation->solid_faces_end(); ++fit) {

        FaceHandle face = fit;
        VectorS2 position = dual_vertex(face);

        std::vector<Arc> vertex_arcs;
        for (int i = 0; i < 3; ++i) {
            Arc arc = dual_arc(Edge(face, i));

            if ((arc.target() - position).squaredNorm() < (arc.source() - position).squaredNorm()) {
                arc = Arc(arc.target(), arc.source(), -arc.normal());
            }

            vertex_arcs.push_back(arc);
        }

        VectorS2 radial = position.normalized();
        VectorS2 ref = (std::abs(radial.z()) < 0.9)
            ? VectorS2(0, 0, 1).cross(radial).normalized()
            : VectorS2(1, 0, 0).cross(radial).normalized();
        VectorS2 perp = radial.cross(ref).normalized();

        std::sort(vertex_arcs.begin(), vertex_arcs.end(),
            [&](const Arc& a, const Arc& b) {
                VectorS2 dir_a = (a.interpolate(0.01) - position).normalized();
                VectorS2 dir_b = (b.interpolate(0.01) - position).normalized();

                double angle_a = std::atan2(dir_a.dot(perp), dir_a.dot(ref));
                double angle_b = std::atan2(dir_b.dot(perp), dir_b.dot(ref));
                return angle_a < angle_b;
            });

        result.push_back({position, std::move(vertex_arcs)});
    }

    return result;
}

inline std::vector<Arc> Sphere::unique_arcs() const {
    if (_triangulation->dimension() < 2) {
        return {};
    }

    std::vector<Arc> result;
    std::set<std::pair<size_t, size_t>> seen_edges;

    for (size_t cell_index = 0; cell_index < size(); ++cell_index) {
        for (const auto& edge_info : cell_edges(cell_index)) {
            size_t a = cell_index;
            size_t b = edge_info.neighbor_index;

            auto edge_key = (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
            if (seen_edges.insert(edge_key).second) {
                result.push_back(edge_info.arc);
            }
        }
    }

    return result;
}

inline size_t Sphere::vertex_index(VertexHandle handle) const {
    auto it = _handle_to_index.find(handle);
    if (it != _handle_to_index.end()) {
        return it->second;
    }
    return size();
}

inline std::vector<CellEdgeInfo> Sphere::cell_edges(size_t index) const {
    if (_triangulation->dimension() < 2) {
        return {};
    }

    std::vector<CellEdgeInfo> result;
    VertexHandle vertex_handle = _handles[index];

    for (const auto& edge : incident_edges_range(vertex_handle)) {
        auto face = edge.first;
        int edge_index = edge.second;

        VertexHandle v1 = face->vertex((edge_index + 1) % 3);
        VertexHandle v2 = face->vertex((edge_index + 2) % 3);

        VertexHandle neighbor_handle = (v1 == vertex_handle) ? v2 : v1;
        size_t neighbor_index = vertex_index(neighbor_handle);

        if (neighbor_index < size()) {
            FaceHandle target_face = face->neighbor(edge_index);

            result.push_back({
                neighbor_index,
                vertex_index(face->vertex(edge_index)),
                vertex_index(target_face->vertex(target_face->index(face))),
                dual_arc(edge)
            });
        }
    }

    return result;
}

inline VectorS2 Sphere::dual_vertex(FaceHandle face) const {
    Vector3 p = to_vector3(_triangulation->point(face->vertex(0)));
    Vector3 q = to_vector3(_triangulation->point(face->vertex(1)));
    Vector3 r = to_vector3(_triangulation->point(face->vertex(2)));
    return (q - p).cross(r - p).normalized();
}

// A Voronoi edge lies on the bisector of its two sites, whose normal is
// known exactly; deriving it from the two Voronoi vertices instead loses
// accuracy in proportion to how short the edge is.
inline Arc Sphere::dual_arc(const Edge& edge) const {
    VectorS2 source = dual_vertex(edge.first);
    VectorS2 target = dual_vertex(edge.first->neighbor(edge.second));
    Vector3 first_site = to_vector3(_triangulation->point(edge.first->vertex((edge.second + 1) % 3)));
    Vector3 second_site = to_vector3(_triangulation->point(edge.first->vertex((edge.second + 2) % 3)));
    VectorS2 bisector_normal = VectorS2(first_site - second_site).normalized();
    Arc arc(source, target, bisector_normal);

    return arc.length() <= M_PI ? arc : Arc(source, target, -bisector_normal);
}

}

#endif //GEOMETRY_ART_VORONOI_SPHERICAL_CORE_SPHERE_HPP_
