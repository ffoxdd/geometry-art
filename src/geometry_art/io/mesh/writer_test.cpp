#include "writer.hpp"
#include "types.hpp"
#include "../../cgal/types.hpp"
#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

namespace geometry_art::io::mesh {
namespace {

SurfaceMesh tetrahedron() {
    SurfaceMesh mesh;
    VertexIndex a = mesh.add_vertex(cgal::Point3(0, 0, 0));
    VertexIndex b = mesh.add_vertex(cgal::Point3(1, 0, 0));
    VertexIndex c = mesh.add_vertex(cgal::Point3(0, 1, 0));
    VertexIndex d = mesh.add_vertex(cgal::Point3(0, 0, 1));

    mesh.add_face(a, c, b);
    mesh.add_face(a, b, d);
    mesh.add_face(b, c, d);
    mesh.add_face(a, d, c);

    return mesh;
}

std::filesystem::path scratch(const std::string& name) {
    std::filesystem::path directory = std::filesystem::temp_directory_path() / "geometry-art-writer-test";
    std::filesystem::create_directories(directory);
    return directory / name;
}

TEST(WriterTest, EveryFormatIsNamedByItsExtension) {
    for (const std::string& name : format_names()) {
        EXPECT_EQ(extension_of(format_named(name)), "." + name);
    }
}

TEST(WriterTest, AnUnknownFormatIsRefused) {
    SurfaceMesh mesh = tetrahedron();

    EXPECT_THROW(format_named("gltf"), std::invalid_argument);
    EXPECT_THROW(write(mesh, scratch("model.gltf")), std::invalid_argument);
}

TEST(WriterTest, EveryFormatWritesAFile) {
    for (const std::string& name : format_names()) {
        SurfaceMesh mesh = tetrahedron();
        std::filesystem::path path = scratch("model." + name);
        std::filesystem::remove(path);

        write(mesh, path);

        EXPECT_GT(std::filesystem::file_size(path), 0u) << name;
    }
}

TEST(WriterTest, StlIsBinaryWithOneRecordPerFacet) {
    SurfaceMesh mesh = tetrahedron();
    std::filesystem::path path = scratch("model.stl");

    write(mesh, path);

    EXPECT_EQ(std::filesystem::file_size(path), 84u + 50u * mesh.number_of_faces());
}

} // namespace
} // namespace geometry_art::io::mesh
