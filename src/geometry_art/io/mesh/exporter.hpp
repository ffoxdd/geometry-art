#ifndef GEOMETRY_ART_IO_MESH_EXPORTER_HPP_
#define GEOMETRY_ART_IO_MESH_EXPORTER_HPP_

#include "voronoi_sphere_wireframe_builder.hpp"
#include "../../voronoi/spherical/core/sphere.hpp"
#include <CGAL/IO/PLY.h>
#include <CGAL/IO/STL.h>
#include <fstream>
#include <stdexcept>
#include <format>
#include <string>

namespace geometry_art::io::mesh {

using voronoi::spherical::Sphere;

enum class Format { PLY, STL };

class Exporter {
public:
    static void save(
        const Sphere &sphere,
        const std::string &filename,
        Format format = Format::STL,
        double arc_thickness = 0.02,
        double max_edge_length = 0.05,
        double scale = 1.0
    );
};

inline void Exporter::save(
    const Sphere &sphere,
    const std::string &filename,
    Format format,
    double arc_thickness,
    double max_edge_length,
    double scale
) {
    VoronoiSphereWireframeBuilder builder(arc_thickness, max_edge_length);
    SurfaceMesh mesh = builder.build(sphere);

    if (scale != 1.0) {
        for (auto vertex : mesh.vertices()) {
            auto& point = mesh.point(vertex);
            point = cgal::Point3(
                point.x() * scale,
                point.y() * scale,
                point.z() * scale
            );
        }
    }

    std::ofstream stream(filename);

    if (!stream) {
        throw std::runtime_error(std::format("Cannot open file for writing: {}", filename));
    }

    bool success = false;

    switch (format) {
        case Format::PLY:
            success = CGAL::IO::write_PLY(stream, mesh);
            break;
        case Format::STL:
            success = CGAL::IO::write_STL(stream, mesh);
            break;
    }

    if (!success) {
        throw std::runtime_error(std::format("Cannot write file: {}", filename));
    }
}

} // namespace geometry_art::io::mesh

#endif //GEOMETRY_ART_IO_MESH_EXPORTER_HPP_
