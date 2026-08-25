#include "gtest/gtest.h"
#include "ranges.hpp"
#include <vector>
#include <list>
#include <array>

using namespace geometry_art::std_ext;

TEST(CircularAdjacentPairsTest, EmptyRange) {
    std::vector<int> empty;
    auto pairs = circular_adjacent_pairs(empty);

    EXPECT_EQ(pairs.begin(), pairs.end());
}

TEST(CircularAdjacentPairsTest, SingleElement) {
    std::vector<int> vec{42};
    auto pairs = circular_adjacent_pairs(vec);

    auto it = pairs.begin();
    ASSERT_NE(it, pairs.end());

    auto [first, second] = *it;
    EXPECT_EQ(first, 42);
    EXPECT_EQ(second, 42);

    ++it;
    EXPECT_EQ(it, pairs.end());
}

TEST(CircularAdjacentPairsTest, TwoElements) {
    std::vector<int> vec{1, 2};
    auto pairs = circular_adjacent_pairs(vec);

    auto it = pairs.begin();
    ASSERT_NE(it, pairs.end());

    auto [first, second] = *it;
    EXPECT_EQ(first, 1);
    EXPECT_EQ(second, 2);

    ++it;
    ASSERT_NE(it, pairs.end());

    auto [first2, second2] = *it;
    EXPECT_EQ(first2, 2);
    EXPECT_EQ(second2, 1);

    ++it;
    EXPECT_EQ(it, pairs.end());
}

TEST(CircularAdjacentPairsTest, ThreeElements) {
    std::vector<int> vec{10, 20, 30};
    auto pairs = circular_adjacent_pairs(vec);

    auto it = pairs.begin();

    auto [a1, a2] = *it;
    EXPECT_EQ(a1, 10);
    EXPECT_EQ(a2, 20);
    ++it;

    auto [b1, b2] = *it;
    EXPECT_EQ(b1, 20);
    EXPECT_EQ(b2, 30);
    ++it;

    auto [c1, c2] = *it;
    EXPECT_EQ(c1, 30);
    EXPECT_EQ(c2, 10);
    ++it;

    EXPECT_EQ(it, pairs.end());
}

TEST(CircularAdjacentPairsTest, MultipleElements) {
    std::vector<int> vec{1, 2, 3, 4, 5};
    auto pairs = circular_adjacent_pairs(vec);

    std::vector<std::pair<int, int>> expected{
        {1, 2},
        {2, 3},
        {3, 4},
        {4, 5},
        {5, 1}
    };

    size_t index = 0;
    for (const auto &pair : pairs) {
        ASSERT_LT(index, expected.size());
        EXPECT_EQ(pair.first, expected[index].first);
        EXPECT_EQ(pair.second, expected[index].second);
        index++;
    }

    EXPECT_EQ(index, expected.size());
}

TEST(CircularAdjacentPairsTest, WorksWithList) {
    std::list<std::string> list{"a", "b", "c"};
    auto pairs = circular_adjacent_pairs(list);

    auto it = pairs.begin();

    auto [a1, a2] = *it;
    EXPECT_EQ(a1, "a");
    EXPECT_EQ(a2, "b");
    ++it;

    auto [b1, b2] = *it;
    EXPECT_EQ(b1, "b");
    EXPECT_EQ(b2, "c");
    ++it;

    auto [c1, c2] = *it;
    EXPECT_EQ(c1, "c");
    EXPECT_EQ(c2, "a");
    ++it;

    EXPECT_EQ(it, pairs.end());
}

TEST(CircularAdjacentPairsTest, WorksWithArray) {
    std::array<int, 4> arr{100, 200, 300, 400};
    auto pairs = circular_adjacent_pairs(arr);

    std::vector<std::pair<int, int>> expected{
        {100, 200},
        {200, 300},
        {300, 400},
        {400, 100}
    };

    size_t index = 0;
    for (const auto &pair : pairs) {
        ASSERT_LT(index, expected.size());
        EXPECT_EQ(pair.first, expected[index].first);
        EXPECT_EQ(pair.second, expected[index].second);
        index++;
    }

    EXPECT_EQ(index, expected.size());
}


TEST(CircularAdjacentPairsTest, WorksWithConstRange) {
    const std::vector<int> vec{5, 10, 15};
    auto pairs = circular_adjacent_pairs(vec);

    auto it = pairs.begin();
    auto [a1, a2] = *it;
    EXPECT_EQ(a1, 5);
    EXPECT_EQ(a2, 10);
}

