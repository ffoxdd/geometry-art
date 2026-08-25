#include "gtest/gtest.h"
#include "indexed_kd_tree.hpp"
#include <boost/iterator/counting_iterator.hpp>

using namespace geometry_art::geometry::spherical;
namespace cgal = geometry_art::cgal;

TEST(IndexedPointMapTest, ReturnsCorrectPointForIndex) {
    std::vector<cgal::Point3> points = {
        cgal::Point3(1, 0, 0),
        cgal::Point3(0, 1, 0),
        cgal::Point3(0, 0, 1)
    };

    IndexedPointMap<std::vector<cgal::Point3>> map(points);

    EXPECT_EQ(map[0], points[0]);
    EXPECT_EQ(map[1], points[1]);
    EXPECT_EQ(map[2], points[2]);
}

TEST(IndexedPointMapTest, GetFunctionReturnsCorrectPoint) {
    std::vector<cgal::Point3> points = {
        cgal::Point3(1, 2, 3),
        cgal::Point3(4, 5, 6)
    };

    IndexedPointMap<std::vector<cgal::Point3>> map(points);

    EXPECT_EQ(get(map, 0), points[0]);
    EXPECT_EQ(get(map, 1), points[1]);
}

TEST(IndexedKDTreeTest, QueriesReturnIndices) {
    std::vector<cgal::Point3> points = {
        cgal::Point3(1, 0, 0),
        cgal::Point3(0, 1, 0),
        cgal::Point3(0, 0, 1),
        cgal::Point3(-1, 0, 0)
    };

    IndexedPointMap<std::vector<cgal::Point3>> point_map(points);
    IndexedSearchTraits search_traits(point_map);

    IndexedKDTree tree(
        boost::counting_iterator<std::size_t>(0),
        boost::counting_iterator<std::size_t>(points.size()),
        IndexedKDTree::Splitter(),
        search_traits
    );
    tree.build();

    // Query near point 0 (1,0,0) with small radius
    IndexedFuzzySphere query(points[0], 0.1, 0.0, search_traits);
    std::vector<size_t> results;
    tree.search(std::back_inserter(results), query);

    ASSERT_EQ(results.size(), 1);
    EXPECT_EQ(results[0], 0);
}

TEST(IndexedKDTreeTest, QueryFindsMultipleNeighbors) {
    std::vector<cgal::Point3> points = {
        cgal::Point3(1, 0, 0),
        cgal::Point3(0.9, 0.1, 0),
        cgal::Point3(0.9, -0.1, 0),
        cgal::Point3(-1, 0, 0)
    };

    IndexedPointMap<std::vector<cgal::Point3>> point_map(points);
    IndexedSearchTraits search_traits(point_map);

    IndexedKDTree tree(
        boost::counting_iterator<std::size_t>(0),
        boost::counting_iterator<std::size_t>(points.size()),
        IndexedKDTree::Splitter(),
        search_traits
    );
    tree.build();

    // Query centered at point 0 with radius that includes points 1 and 2
    IndexedFuzzySphere query(points[0], 0.25, 0.0, search_traits);
    std::vector<size_t> results;
    tree.search(std::back_inserter(results), query);

    EXPECT_EQ(results.size(), 3);

    std::sort(results.begin(), results.end());
    EXPECT_EQ(results[0], 0);
    EXPECT_EQ(results[1], 1);
    EXPECT_EQ(results[2], 2);
}
