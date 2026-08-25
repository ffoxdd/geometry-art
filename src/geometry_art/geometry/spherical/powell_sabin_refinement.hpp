#ifndef GEOMETRY_ART_GEOMETRY_SPHERICAL_POWELL_SABIN_REFINEMENT_HPP_
#define GEOMETRY_ART_GEOMETRY_SPHERICAL_POWELL_SABIN_REFINEMENT_HPP_

#include "triangle_mesh.hpp"
#include "../../types.hpp"
#include <CGAL/assertions.h>
#include <array>
#include <cstddef>
#include <map>
#include <utility>
#include <vector>

namespace geometry_art::geometry::spherical {

// Splits every triangle of a mesh into six, the refinement a C1 quadratic
// spline needs.
//
// The split point on a shared edge is where the great circle joining the two
// neighbouring interior points crosses it. That placement is the whole trick:
// it makes the two sides of the edge agree on the direction the surface
// leaves it in, which is what allows one quadratic per piece to be smooth
// across the seam. Any other split point gives a mesh on which no such
// spline exists.
class PowellSabinRefinement {
 public:
    enum Point {
        CORNER_0 = 0,
        CORNER_1 = 1,
        CORNER_2 = 2,
        SPLIT_01 = 3,
        SPLIT_12 = 4,
        SPLIT_20 = 5,
        INTERIOR = 6,
        POINT_COUNT = 7
    };

    struct Cell {
        std::array<VectorS2, POINT_COUNT> points;
        std::array<size_t, 3> corner_vertices;
        std::array<size_t, POINT_COUNT> refined_vertices;
    };

    // Sub-triangles as indices into a cell's points, walking the boundary so
    // that each keeps the orientation of the triangle it came from.
    static constexpr std::array<std::array<int, 3>, 6> SUB_TRIANGLES = {{
        {CORNER_0, SPLIT_01, INTERIOR},
        {SPLIT_01, CORNER_1, INTERIOR},
        {CORNER_1, SPLIT_12, INTERIOR},
        {SPLIT_12, CORNER_2, INTERIOR},
        {CORNER_2, SPLIT_20, INTERIOR},
        {SPLIT_20, CORNER_0, INTERIOR}
    }};

    // Pairs of sub-triangles sharing an interior edge, where the smoothness
    // conditions inside a cell are written.
    static constexpr std::array<std::array<int, 2>, 6> ADJACENT_SUB_TRIANGLES = {{
        {5, 0}, {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}
    }};

    explicit PowellSabinRefinement(const TriangleMesh& mesh);

    [[nodiscard]] const std::vector<Cell>& cells() const { return _cells; }
    [[nodiscard]] const TriangleMesh& mesh() const { return _mesh; }

 private:
    using Edge = std::pair<size_t, size_t>;

    std::vector<Cell> _cells;
    TriangleMesh _mesh;

    [[nodiscard]] static std::vector<VectorS2> interior_points(const TriangleMesh& mesh);
    [[nodiscard]] static std::map<Edge, std::vector<size_t>> edge_triangles(const TriangleMesh& mesh);
    [[nodiscard]] static VectorS2 split_point(
        const VectorS2& from,
        const VectorS2& to,
        const std::vector<size_t>& triangles,
        const std::vector<VectorS2>& interiors
    );
    [[nodiscard]] static Edge edge_of(size_t first, size_t second);
};

inline PowellSabinRefinement::PowellSabinRefinement(const TriangleMesh& mesh) {
    std::vector<VectorS2> interiors = interior_points(mesh);
    std::map<Edge, std::vector<size_t>> incident = edge_triangles(mesh);

    _mesh.vertices = mesh.vertices;

    std::map<Edge, size_t> split_vertices;

    for (const auto& [edge, triangles] : incident) {
        split_vertices[edge] = _mesh.vertices.size();
        _mesh.vertices.push_back(
            split_point(mesh.vertices[edge.first], mesh.vertices[edge.second], triangles, interiors)
        );
    }

    _cells.reserve(mesh.triangles.size());

    for (size_t index = 0; index < mesh.triangles.size(); ++index) {
        const std::array<size_t, 3>& corners = mesh.triangles[index];

        Cell cell;
        cell.corner_vertices = corners;
        cell.refined_vertices[CORNER_0] = corners[0];
        cell.refined_vertices[CORNER_1] = corners[1];
        cell.refined_vertices[CORNER_2] = corners[2];
        cell.refined_vertices[SPLIT_01] = split_vertices.at(edge_of(corners[0], corners[1]));
        cell.refined_vertices[SPLIT_12] = split_vertices.at(edge_of(corners[1], corners[2]));
        cell.refined_vertices[SPLIT_20] = split_vertices.at(edge_of(corners[2], corners[0]));

        cell.refined_vertices[INTERIOR] = _mesh.vertices.size();
        _mesh.vertices.push_back(interiors[index]);

        for (int point = 0; point < POINT_COUNT; ++point) {
            cell.points[point] = _mesh.vertices[cell.refined_vertices[point]];
        }

        for (const std::array<int, 3>& sub_triangle : SUB_TRIANGLES) {
            _mesh.triangles.push_back({
                cell.refined_vertices[sub_triangle[0]],
                cell.refined_vertices[sub_triangle[1]],
                cell.refined_vertices[sub_triangle[2]]
            });
        }

        _cells.push_back(cell);
    }
}

inline std::vector<VectorS2> PowellSabinRefinement::interior_points(const TriangleMesh& mesh) {
    std::vector<VectorS2> result;
    result.reserve(mesh.triangles.size());

    for (const std::array<size_t, 3>& triangle : mesh.triangles) {
        result.push_back(
            VectorS2(mesh.vertices[triangle[0]] + mesh.vertices[triangle[1]] + mesh.vertices[triangle[2]]).normalized()
        );
    }

    return result;
}

inline std::map<PowellSabinRefinement::Edge, std::vector<size_t>>
PowellSabinRefinement::edge_triangles(const TriangleMesh& mesh) {
    std::map<Edge, std::vector<size_t>> result;

    for (size_t index = 0; index < mesh.triangles.size(); ++index) {
        const std::array<size_t, 3>& triangle = mesh.triangles[index];

        for (int corner = 0; corner < 3; ++corner) {
            result[edge_of(triangle[corner], triangle[(corner + 1) % 3])].push_back(index);
        }
    }

    return result;
}

inline VectorS2 PowellSabinRefinement::split_point(
    const VectorS2& from,
    const VectorS2& to,
    const std::vector<size_t>& triangles,
    const std::vector<VectorS2>& interiors
) {
    VectorS2 midpoint = VectorS2(from + to).normalized();

    if (triangles.size() != 2) {
        return midpoint;
    }

    Vector3 through_interiors = interiors[triangles[0]].cross(interiors[triangles[1]]);
    Vector3 along_edge = from.cross(to);
    Vector3 crossing = through_interiors.cross(along_edge);

    if (crossing.norm() < GEOMETRIC_EPSILON) {
        return midpoint;
    }

    VectorS2 candidate = crossing.normalized();

    if (candidate.dot(midpoint) < 0.0) {
        candidate = -candidate;
    }

    // Outside the edge there is no spline to build, so a mesh that would
    // produce one is a precondition failure rather than something to paper
    // over with the midpoint.
    CGAL_postcondition(candidate.dot(from) > 0.0 && candidate.dot(to) > 0.0);

    return candidate;
}

inline PowellSabinRefinement::Edge PowellSabinRefinement::edge_of(size_t first, size_t second) {
    return first < second ? Edge(first, second) : Edge(second, first);
}

} // namespace geometry_art::geometry::spherical

#endif //GEOMETRY_ART_GEOMETRY_SPHERICAL_POWELL_SABIN_REFINEMENT_HPP_
