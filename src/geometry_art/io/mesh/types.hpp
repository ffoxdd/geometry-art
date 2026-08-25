#ifndef GEOMETRY_ART_IO_MESH_TYPES_HPP_
#define GEOMETRY_ART_IO_MESH_TYPES_HPP_

#include "../../cgal/types.hpp"
#include <CGAL/Surface_mesh.h>

namespace geometry_art::io::mesh {

using SurfaceMesh = ::CGAL::Surface_mesh<cgal::Point3>;
using VertexIndex = SurfaceMesh::Vertex_index;

} // namespace geometry_art::io::mesh

#endif //GEOMETRY_ART_IO_MESH_TYPES_HPP_
