#ifndef GLOBEART_SRC_GLOBE_IO_MESH_VORONOI_SPHERE_WIREFRAME_BUILDER_HPP_
#define GLOBEART_SRC_GLOBE_IO_MESH_VORONOI_SPHERE_WIREFRAME_BUILDER_HPP_

#include "types.hpp"
#include "../../cgal/types.hpp"
#include "../../types.hpp"
#include "../../geometry/spherical/arc.hpp"
#include "../../voronoi/spherical/core/sphere.hpp"
#include <map>
#include <set>
#include <vector>
#include <cmath>

namespace globe::io::mesh {

using geometry::spherical::Arc;
using voronoi::spherical::Sphere;
using voronoi::spherical::VoronoiVertex;

class VoronoiSphereWireframeBuilder {
public:
    VoronoiSphereWireframeBuilder(double thickness, double max_edge_length = 0.05);

    SurfaceMesh build(const Sphere& sphere);

private:
    struct CellVertex {
        VectorS2 position;
        VectorS2 chamfer;
    };

    struct CellData {
        std::vector<CellVertex> vertices;
    };

    void compute_all_cells(const Sphere& sphere);
    void build_all_bars(const Sphere& sphere);
    void build_vertex_caps(const Sphere& sphere);

    VectorS2 compute_chamfer_point(
        const VectorS2& prev,
        const VectorS2& curr,
        const VectorS2& next
    ) const;

    VectorS2 slerp(const VectorS2& a, const VectorS2& b, double t) const;
    double arc_length(const VectorS2& a, const VectorS2& b) const;
    int compute_sample_count(
        const VectorS2& a_start, const VectorS2& a_end,
        const VectorS2& b_start, const VectorS2& b_end
    ) const;

    void build_bar_complete(
        const VectorS2& chamfer_a_start,
        const VectorS2& chamfer_a_end,
        const VectorS2& chamfer_b_start,
        const VectorS2& chamfer_b_end
    );

    void build_flat_cap(
        const std::vector<VectorS2>& polygon_vertices,
        bool outward_normal
    );

    VertexIndex get_or_create_vertex(const VectorS2& point);

    double _thickness;
    double _half_thickness;
    double _max_edge_length;
    SurfaceMesh _mesh;

    std::map<std::tuple<double, double, double>, VertexIndex> _vertex_cache;
    std::vector<CellData> _cells;
};

// Implementation

inline VoronoiSphereWireframeBuilder::VoronoiSphereWireframeBuilder(
    double thickness,
    double max_edge_length
) : _thickness(thickness),
    _half_thickness(thickness / 2.0),
    _max_edge_length(max_edge_length) {
}

inline SurfaceMesh VoronoiSphereWireframeBuilder::build(const Sphere& sphere) {
    _mesh = SurfaceMesh();
    _vertex_cache.clear();
    _cells.clear();

    // Need at least 3 points for a valid Voronoi diagram on a sphere
    if (sphere.size() < 3) {
        return SurfaceMesh();
    }

    compute_all_cells(sphere);
    build_all_bars(sphere);
    build_vertex_caps(sphere);

    SurfaceMesh result = std::move(_mesh);
    _mesh = SurfaceMesh();
    _vertex_cache.clear();
    _cells.clear();

    return result;
}

inline void VoronoiSphereWireframeBuilder::compute_all_cells(const Sphere& sphere) {
    _cells.resize(sphere.size());

    for (size_t cell_index = 0; cell_index < sphere.size(); ++cell_index) {
        auto polygon = sphere.cell(cell_index);
        const auto& arcs = polygon.arcs();
        size_t n = arcs.size();
        if (n < 3) continue;

        CellData& cell_data = _cells[cell_index];

        for (size_t i = 0; i < n; ++i) {
            VectorS2 prev = arcs[(i + n - 1) % n].source();
            VectorS2 curr = arcs[i].source();
            VectorS2 next = arcs[(i + 1) % n].source();

            VectorS2 chamfer = compute_chamfer_point(prev, curr, next);
            cell_data.vertices.push_back({curr, chamfer});
        }
    }
}

inline void VoronoiSphereWireframeBuilder::build_all_bars(const Sphere& sphere) {
    std::set<std::pair<size_t, size_t>> processed_edges;

    for (size_t cell_index = 0; cell_index < sphere.size(); ++cell_index) {
        const auto& cell_data = _cells[cell_index];
        size_t n = cell_data.vertices.size();
        if (n < 3) continue;

        auto cell_edges = sphere.cell_edges(cell_index);

        for (size_t i = 0; i < n; ++i) {
            size_t next_i = (i + 1) % n;
            VectorS2 vertex_start = cell_data.vertices[i].position;
            VectorS2 vertex_end = cell_data.vertices[next_i].position;

            size_t neighbor_index = SIZE_MAX;
            for (const auto& edge_info : cell_edges) {
                VectorS2 arc_start = edge_info.arc.source();
                VectorS2 arc_end = edge_info.arc.target();

                bool match_forward = (vertex_start - arc_start).norm() < 1e-6 &&
                                     (vertex_end - arc_end).norm() < 1e-6;
                bool match_reverse = (vertex_start - arc_end).norm() < 1e-6 &&
                                     (vertex_end - arc_start).norm() < 1e-6;

                if (match_forward || match_reverse) {
                    neighbor_index = edge_info.neighbor_index;
                    break;
                }
            }

            if (neighbor_index == SIZE_MAX) continue;

            auto edge_key = (cell_index < neighbor_index)
                ? std::make_pair(cell_index, neighbor_index)
                : std::make_pair(neighbor_index, cell_index);

            if (processed_edges.count(edge_key)) continue;
            processed_edges.insert(edge_key);

            VectorS2 chamfer_a_start = cell_data.vertices[i].chamfer;
            VectorS2 chamfer_a_end = cell_data.vertices[next_i].chamfer;

            const auto& neighbor_data = _cells[neighbor_index];
            VectorS2 chamfer_b_start, chamfer_b_end;

            for (size_t j = 0; j < neighbor_data.vertices.size(); ++j) {
                if ((neighbor_data.vertices[j].position - vertex_start).norm() < 1e-6) {
                    chamfer_b_start = neighbor_data.vertices[j].chamfer;
                }
                if ((neighbor_data.vertices[j].position - vertex_end).norm() < 1e-6) {
                    chamfer_b_end = neighbor_data.vertices[j].chamfer;
                }
            }

            build_bar_complete(chamfer_a_start, chamfer_a_end, chamfer_b_start, chamfer_b_end);
        }
    }
}

inline void VoronoiSphereWireframeBuilder::build_vertex_caps(const Sphere& sphere) {
    auto voronoi_vertices = sphere.vertices();
    double outer_radius = 1.0 + _half_thickness;
    double inner_radius = 1.0 - _half_thickness;

    for (const auto& vertex : voronoi_vertices) {
        VectorS2 vertex_pos = vertex.position;

        std::vector<VectorS2> chamfer_points;
        for (size_t cell_index = 0; cell_index < _cells.size(); ++cell_index) {
            const auto& cell_data = _cells[cell_index];
            for (const auto& cv : cell_data.vertices) {
                if ((cv.position - vertex_pos).norm() < 1e-6) {
                    chamfer_points.push_back(cv.chamfer);
                    break;
                }
            }
        }

        if (chamfer_points.size() < 3) continue;

        VectorS2 radial = vertex_pos.normalized();
        VectorS2 ref = (std::abs(radial.z()) < 0.9)
            ? VectorS2(0, 0, 1).cross(radial).normalized()
            : VectorS2(1, 0, 0).cross(radial).normalized();

        std::sort(chamfer_points.begin(), chamfer_points.end(),
            [&](const VectorS2& a, const VectorS2& b) {
                VectorS2 da = (a - vertex_pos).normalized();
                VectorS2 db = (b - vertex_pos).normalized();
                double angle_a = std::atan2(radial.cross(ref).dot(da), ref.dot(da));
                double angle_b = std::atan2(radial.cross(ref).dot(db), ref.dot(db));
                return angle_a < angle_b;
            });

        std::vector<VectorS2> outer_points, inner_points;
        for (const auto& chamfer : chamfer_points) {
            outer_points.push_back(chamfer * outer_radius);
            inner_points.push_back(chamfer * inner_radius);
        }

        build_flat_cap(outer_points, true);
        build_flat_cap(inner_points, false);
    }
}

inline VectorS2 VoronoiSphereWireframeBuilder::compute_chamfer_point(
    const VectorS2& prev,
    const VectorS2& curr,
    const VectorS2& next
) const {
    VectorS2 to_prev = (prev - curr).normalized();
    VectorS2 to_next = (next - curr).normalized();

    VectorS2 bisector = to_prev + to_next;
    double bisector_len = bisector.norm();

    if (bisector_len < 1e-10) {
        VectorS2 radial = curr.normalized();
        bisector = radial.cross(to_next).normalized();
    } else {
        bisector = bisector / bisector_len;
    }

    double cos_angle = to_prev.dot(to_next);
    double angle = std::acos(std::clamp(cos_angle, -1.0, 1.0));
    double half_angle = angle / 2.0;
    double sin_half = std::sin(half_angle);

    double distance = _half_thickness / std::max(sin_half, 0.1);

    VectorS2 chamfer = curr + bisector * distance;
    return chamfer.normalized();
}

inline VectorS2 VoronoiSphereWireframeBuilder::slerp(
    const VectorS2& a,
    const VectorS2& b,
    double t
) const {
    double dot = a.dot(b);
    dot = std::clamp(dot, -1.0, 1.0);
    double angle = std::acos(dot);
    if (angle < 1e-10) return a;
    double sin_angle = std::sin(angle);
    return a * (std::sin((1.0 - t) * angle) / sin_angle) +
           b * (std::sin(t * angle) / sin_angle);
}

inline double VoronoiSphereWireframeBuilder::arc_length(
    const VectorS2& a,
    const VectorS2& b
) const {
    double dot = a.normalized().dot(b.normalized());
    dot = std::clamp(dot, -1.0, 1.0);
    return std::acos(dot);
}

inline int VoronoiSphereWireframeBuilder::compute_sample_count(
    const VectorS2& a_start, const VectorS2& a_end,
    const VectorS2& b_start, const VectorS2& b_end
) const {
    double length_a = arc_length(a_start, a_end);
    double length_b = arc_length(b_start, b_end);
    double max_length = std::max(length_a, length_b);

    int samples = static_cast<int>(std::ceil(max_length / _max_edge_length));
    return std::max(samples, 1);
}

inline void VoronoiSphereWireframeBuilder::build_bar_complete(
    const VectorS2& chamfer_a_start,
    const VectorS2& chamfer_a_end,
    const VectorS2& chamfer_b_start,
    const VectorS2& chamfer_b_end
) {
    int sample_count = compute_sample_count(
        chamfer_a_start, chamfer_a_end,
        chamfer_b_start, chamfer_b_end
    );

    for (int s = 0; s < sample_count; ++s) {
        double t0 = static_cast<double>(s) / sample_count;
        double t1 = static_cast<double>(s + 1) / sample_count;

        VectorS2 point_a0 = slerp(chamfer_a_start, chamfer_a_end, t0);
        VectorS2 point_a1 = slerp(chamfer_a_start, chamfer_a_end, t1);
        VectorS2 point_b0 = slerp(chamfer_b_start, chamfer_b_end, t0);
        VectorS2 point_b1 = slerp(chamfer_b_start, chamfer_b_end, t1);

        VertexIndex vertex_a0_outer = get_or_create_vertex(point_a0 * (1.0 + _half_thickness));
        VertexIndex vertex_a0_inner = get_or_create_vertex(point_a0 * (1.0 - _half_thickness));
        VertexIndex vertex_a1_outer = get_or_create_vertex(point_a1 * (1.0 + _half_thickness));
        VertexIndex vertex_a1_inner = get_or_create_vertex(point_a1 * (1.0 - _half_thickness));
        VertexIndex vertex_b0_outer = get_or_create_vertex(point_b0 * (1.0 + _half_thickness));
        VertexIndex vertex_b0_inner = get_or_create_vertex(point_b0 * (1.0 - _half_thickness));
        VertexIndex vertex_b1_outer = get_or_create_vertex(point_b1 * (1.0 + _half_thickness));
        VertexIndex vertex_b1_inner = get_or_create_vertex(point_b1 * (1.0 - _half_thickness));

        // Top face (outer surface)
        _mesh.add_face(vertex_a0_outer, vertex_b0_outer, vertex_b1_outer);
        _mesh.add_face(vertex_a0_outer, vertex_b1_outer, vertex_a1_outer);

        // Bottom face (inner surface)
        _mesh.add_face(vertex_a0_inner, vertex_a1_inner, vertex_b1_inner);
        _mesh.add_face(vertex_a0_inner, vertex_b1_inner, vertex_b0_inner);

        // Side A (facing cell A)
        _mesh.add_face(vertex_a0_outer, vertex_a1_outer, vertex_a1_inner);
        _mesh.add_face(vertex_a0_outer, vertex_a1_inner, vertex_a0_inner);

        // Side B (facing cell B)
        _mesh.add_face(vertex_b0_outer, vertex_b0_inner, vertex_b1_inner);
        _mesh.add_face(vertex_b0_outer, vertex_b1_inner, vertex_b1_outer);
    }
}

inline void VoronoiSphereWireframeBuilder::build_flat_cap(
    const std::vector<VectorS2>& polygon_vertices,
    bool outward_normal
) {
    size_t n = polygon_vertices.size();
    if (n < 3) return;

    VectorS2 centroid(0, 0, 0);
    for (const auto& v : polygon_vertices) {
        centroid = centroid + v;
    }
    centroid = centroid * (1.0 / n);

    VertexIndex center_index = get_or_create_vertex(centroid);

    for (size_t i = 0; i < n; ++i) {
        size_t next_i = (i + 1) % n;
        VertexIndex i1 = get_or_create_vertex(polygon_vertices[i]);
        VertexIndex i2 = get_or_create_vertex(polygon_vertices[next_i]);

        if (outward_normal) {
            _mesh.add_face(center_index, i1, i2);
        } else {
            _mesh.add_face(center_index, i2, i1);
        }
    }
}

inline VertexIndex VoronoiSphereWireframeBuilder::get_or_create_vertex(const VectorS2& point) {
    auto key = std::make_tuple(
        std::round(point.x() * 1e6) / 1e6,
        std::round(point.y() * 1e6) / 1e6,
        std::round(point.z() * 1e6) / 1e6
    );

    auto it = _vertex_cache.find(key);
    if (it != _vertex_cache.end()) {
        return it->second;
    }

    cgal::Point3 cgal_point(point.x(), point.y(), point.z());
    VertexIndex index = _mesh.add_vertex(cgal_point);
    _vertex_cache[key] = index;
    return index;
}

} // namespace globe::io::mesh

#endif //GLOBEART_SRC_GLOBE_IO_MESH_VORONOI_SPHERE_WIREFRAME_BUILDER_HPP_
