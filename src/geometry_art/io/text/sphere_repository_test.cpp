#include "sphere_repository.hpp"
#include <gtest/gtest.h>
#include <cmath>
#include <filesystem>

namespace geometry_art::io::text {
namespace {

TEST(SphereRepositoryTest, RoundTrip) {
    Sphere original;
    original.insert(cgal::Point3(1, 0, 0));
    original.insert(cgal::Point3(0, 1, 0));
    original.insert(cgal::Point3(0, 0, 1));
    original.insert(cgal::Point3(-1, 0, 0));
    original.insert(cgal::Point3(0, -1, 0));
    original.insert(cgal::Point3(0, 0, -1));

    std::string path = "/tmp/test_sphere_repository.txt";
    SphereRepository::save(original, path);

    auto loaded = SphereRepository::load(path);

    ASSERT_EQ(loaded->size(), original.size());

    for (size_t i = 0; i < original.size(); ++i) {
        cgal::Point3 orig_site = original.site(i);
        cgal::Point3 load_site = loaded->site(i);

        EXPECT_DOUBLE_EQ(orig_site.x(), load_site.x());
        EXPECT_DOUBLE_EQ(orig_site.y(), load_site.y());
        EXPECT_DOUBLE_EQ(orig_site.z(), load_site.z());
    }

    std::filesystem::remove(path);
}

TEST(SphereRepositoryTest, EmptySphere) {
    Sphere original;

    std::string path = "/tmp/test_sphere_repository_empty.txt";
    SphereRepository::save(original, path);

    auto loaded = SphereRepository::load(path);

    EXPECT_EQ(loaded->size(), 0);

    std::filesystem::remove(path);
}

TEST(SphereRepositoryTest, HighPrecision) {
    Sphere original;
    double inv_sqrt3 = 0.57735026918962576;
    original.insert(cgal::Point3(inv_sqrt3, inv_sqrt3, inv_sqrt3));
    original.insert(cgal::Point3(-inv_sqrt3, -inv_sqrt3, -inv_sqrt3));
    original.insert(cgal::Point3(inv_sqrt3, -inv_sqrt3, -inv_sqrt3));
    original.insert(cgal::Point3(-inv_sqrt3, inv_sqrt3, inv_sqrt3));

    std::string path = "/tmp/test_sphere_repository_precision.txt";
    SphereRepository::save(original, path);

    auto loaded = SphereRepository::load(path);

    ASSERT_EQ(loaded->size(), original.size());

    cgal::Point3 orig_site = original.site(0);
    cgal::Point3 load_site = loaded->site(0);

    EXPECT_DOUBLE_EQ(orig_site.x(), load_site.x());
    EXPECT_DOUBLE_EQ(orig_site.y(), load_site.y());
    EXPECT_DOUBLE_EQ(orig_site.z(), load_site.z());

    std::filesystem::remove(path);
}

TEST(SphereRepositoryTest, FileNotFound) {
    EXPECT_THROW(SphereRepository::load("/nonexistent/path/file.txt"), std::runtime_error);
}

TEST(SphereRepositoryTest, IgnoresComments) {
    std::string path = "/tmp/test_sphere_repository_comments.txt";

    std::ofstream file(path);
    file << "# This is a comment\n";
    file << "# Points: 4\n";
    file << "1 0 0\n";
    file << "# Another comment\n";
    file << "-1 0 0\n";
    file << "0 1 0\n";
    file << "0 -1 0\n";
    file.close();

    auto loaded = SphereRepository::load(path);

    EXPECT_EQ(loaded->size(), 4);

    std::filesystem::remove(path);
}

} // namespace
} // namespace geometry_art::io::text
