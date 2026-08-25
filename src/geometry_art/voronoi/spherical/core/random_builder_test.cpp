#include <gtest/gtest.h>
#include "random_builder.hpp"
#include "../../../testing/mocks/point_generator.hpp"

using namespace geometry_art;
using namespace geometry_art::voronoi;
using namespace geometry_art::voronoi::spherical;
using geometry_art::testing::mocks::MockPointGenerator;

TEST(RandomBuilderTest, CreatesSphereWithCorrectSize) {
    RandomBuilder builder;

    auto sphere = builder.build(10);
    EXPECT_EQ(sphere->size(), 10);
}

TEST(RandomBuilderTest, UsesProvidedGenerator) {
    VectorS2 fixed_point(1, 0, 0);
    MockPointGenerator generator({fixed_point});
    RandomBuilder<MockPointGenerator> builder(generator);

    auto sphere = builder.build(1);

    EXPECT_EQ(sphere->size(), 1);

    cgal::Point3 site = sphere->site(0);
    EXPECT_DOUBLE_EQ(site.x(), fixed_point.x());
    EXPECT_DOUBLE_EQ(site.y(), fixed_point.y());
    EXPECT_DOUBLE_EQ(site.z(), fixed_point.z());
}

