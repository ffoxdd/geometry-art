#ifndef GEOMETRY_ART_IO_MESH_WRITER_HPP_
#define GEOMETRY_ART_IO_MESH_WRITER_HPP_

#include "types.hpp"
#include <CGAL/IO/polygon_mesh_io.h>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

namespace geometry_art::io::mesh {

// The file formats a mesh can be written in, named by their extensions. STL
// is what printers and slicers take; OBJ and PLY keep shared vertices for
// modelling tools; OFF is the plain interchange form.
enum class Format {
    stl,
    obj,
    ply,
    off
};

inline const std::vector<std::string>& format_names() {
    static const std::vector<std::string> names{"stl", "obj", "ply", "off"};
    return names;
}

[[nodiscard]] inline Format format_named(const std::string& name) {
    if (name == "stl") return Format::stl;
    if (name == "obj") return Format::obj;
    if (name == "ply") return Format::ply;
    if (name == "off") return Format::off;

    throw std::invalid_argument("no mesh format is called " + name);
}

[[nodiscard]] inline std::string extension_of(Format format) {
    return "." + format_names()[static_cast<size_t>(format)];
}

// Writes the mesh in the format its path's extension names.
inline void write(SurfaceMesh& mesh, const std::filesystem::path& path) {
    format_named(path.extension().string().substr(path.has_extension() ? 1 : 0));

    if (!CGAL::IO::write_polygon_mesh(path.string(), mesh, CGAL::parameters::stream_precision(10))) {
        throw std::runtime_error("cannot write the mesh to " + path.string());
    }
}

} // namespace geometry_art::io::mesh

#endif //GEOMETRY_ART_IO_MESH_WRITER_HPP_
