#ifndef GLOBEART_SRC_GLOBE_GEOMETRY_SPHERICAL_TRIANGLE_MESH_HPP_
#define GLOBEART_SRC_GLOBE_GEOMETRY_SPHERICAL_TRIANGLE_MESH_HPP_

#include "../../types.hpp"
#include <Eigen/Geometry>
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <utility>
#include <vector>

namespace globe::geometry::spherical {

using globe::VectorS2;

struct TriangleMesh {
    std::vector<VectorS2> vertices;
    std::vector<std::array<size_t, 3>> triangles;

    [[nodiscard]] static TriangleMesh icosphere(int subdivisions);

 private:
    [[nodiscard]] static TriangleMesh icosahedron();
    [[nodiscard]] static TriangleMesh subdivided(const TriangleMesh& mesh);
    static void orient_outward(TriangleMesh& mesh);
};

inline TriangleMesh TriangleMesh::icosphere(int subdivisions) {
    TriangleMesh mesh = icosahedron();

    for (int level = 0; level < subdivisions; ++level) {
        mesh = subdivided(mesh);
    }

    orient_outward(mesh);
    return mesh;
}

inline TriangleMesh TriangleMesh::icosahedron() {
    double phi = (1.0 + std::sqrt(5.0)) / 2.0;
    TriangleMesh mesh;

    for (double a : {-1.0, 1.0}) {
        for (double b : {-phi, phi}) {
            mesh.vertices.push_back(VectorS2(a, b, 0.0).normalized());
            mesh.vertices.push_back(VectorS2(0.0, a, b).normalized());
            mesh.vertices.push_back(VectorS2(b, 0.0, a).normalized());
        }
    }

    mesh.triangles = {
        {0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11},
        {1, 5, 9}, {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
        {3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8}, {3, 8, 9},
        {4, 9, 5}, {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1}
    };

    std::vector<VectorS2> reordered = {
        VectorS2(-1, phi, 0).normalized(), VectorS2(1, phi, 0).normalized(),
        VectorS2(-1, -phi, 0).normalized(), VectorS2(1, -phi, 0).normalized(),
        VectorS2(0, -1, phi).normalized(), VectorS2(0, 1, phi).normalized(),
        VectorS2(0, -1, -phi).normalized(), VectorS2(0, 1, -phi).normalized(),
        VectorS2(phi, 0, -1).normalized(), VectorS2(phi, 0, 1).normalized(),
        VectorS2(-phi, 0, -1).normalized(), VectorS2(-phi, 0, 1).normalized()
    };
    mesh.vertices = reordered;

    return mesh;
}

inline TriangleMesh TriangleMesh::subdivided(const TriangleMesh& mesh) {
    TriangleMesh result;
    result.vertices = mesh.vertices;
    std::map<std::pair<size_t, size_t>, size_t> midpoints;

    auto midpoint = [&](size_t a, size_t b) {
        auto key = a < b ? std::make_pair(a, b) : std::make_pair(b, a);
        auto found = midpoints.find(key);
        if (found != midpoints.end()) {
            return found->second;
        }

        result.vertices.push_back((mesh.vertices[a] + mesh.vertices[b]).normalized());
        size_t index = result.vertices.size() - 1;
        midpoints.emplace(key, index);
        return index;
    };

    for (const auto& [a, b, c] : mesh.triangles) {
        size_t ab = midpoint(a, b);
        size_t bc = midpoint(b, c);
        size_t ca = midpoint(c, a);

        result.triangles.push_back({a, ab, ca});
        result.triangles.push_back({b, bc, ab});
        result.triangles.push_back({c, ca, bc});
        result.triangles.push_back({ab, bc, ca});
    }

    return result;
}

inline void TriangleMesh::orient_outward(TriangleMesh& mesh) {
    for (auto& triangle : mesh.triangles) {
        const VectorS2& a = mesh.vertices[triangle[0]];
        const VectorS2& b = mesh.vertices[triangle[1]];
        const VectorS2& c = mesh.vertices[triangle[2]];

        if ((b - a).cross(c - a).dot(a) < 0.0) {
            std::swap(triangle[1], triangle[2]);
        }
    }
}

} // namespace globe::geometry::spherical

#endif //GLOBEART_SRC_GLOBE_GEOMETRY_SPHERICAL_TRIANGLE_MESH_HPP_
