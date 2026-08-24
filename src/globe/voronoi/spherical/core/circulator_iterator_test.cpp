#include <gtest/gtest.h>
#include "circulator_iterator.hpp"
#include <vector>
#include <algorithm>

using namespace globe::voronoi;
using namespace globe::voronoi::spherical;

namespace {

template<typename T>
class MockCirculator {
 public:
    MockCirculator() : _data(nullptr), _index(0), _size(0) {}

    MockCirculator(std::vector<T>* data, size_t index)
        : _data(data), _index(index), _size(data->size()) {}

    const T& operator*() const {
        return (*_data)[_index];
    }

    const T* operator->() const {
        return &(*_data)[_index];
    }

    MockCirculator& operator++() {
        _index = (_index + 1) % _size;
        return *this;
    }

    friend bool operator==(const MockCirculator& a, const MockCirculator& b) {
        return a._data == b._data && a._index == b._index;
    }

 private:
    std::vector<T>* _data;
    size_t _index;
    size_t _size;
};

}

TEST(CirculatorIteratorTest, IteratesOverAllElements) {
    std::vector<int> data = {1, 2, 3, 4};
    MockCirculator<int> circulator(&data, 0);

    auto begin = CirculatorIterator<MockCirculator<int>, int>(circulator);
    auto end = CirculatorIterator<MockCirculator<int>, int>(circulator);

    std::vector<int> result;
    for (auto it = begin; it != end; ++it) {
        result.push_back(*it);
    }

    EXPECT_EQ(result, data);
}

TEST(CirculatorIteratorTest, DereferenceReturnsCurrentElement) {
    std::vector<int> data = {10, 20, 30};
    MockCirculator<int> circulator(&data, 0);

    CirculatorIterator<MockCirculator<int>, int> it(circulator);

    EXPECT_EQ(*it, 10);
}

TEST(CirculatorIteratorTest, PreIncrementAdvancesIterator) {
    std::vector<int> data = {10, 20, 30};
    MockCirculator<int> circulator(&data, 0);

    CirculatorIterator<MockCirculator<int>, int> it(circulator);

    ++it;

    EXPECT_EQ(*it, 20);
}

TEST(CirculatorIteratorTest, PostIncrementReturnsOldValue) {
    std::vector<int> data = {10, 20, 30};
    MockCirculator<int> circulator(&data, 0);

    CirculatorIterator<MockCirculator<int>, int> it(circulator);

    auto old = it++;

    EXPECT_EQ(*old, 10);
    EXPECT_EQ(*it, 20);
}

TEST(CirculatorIteratorTest, EqualityAtStartReturnsFalse) {
    std::vector<int> data = {1, 2, 3};
    MockCirculator<int> circulator(&data, 0);

    auto a = CirculatorIterator<MockCirculator<int>, int>(circulator);
    auto b = CirculatorIterator<MockCirculator<int>, int>(circulator);

    EXPECT_FALSE(a == b);
}

TEST(CirculatorIteratorTest, EqualityAfterFullCycleReturnsTrue) {
    std::vector<int> data = {1, 2, 3};
    MockCirculator<int> circulator(&data, 0);

    auto begin = CirculatorIterator<MockCirculator<int>, int>(circulator);
    auto end = CirculatorIterator<MockCirculator<int>, int>(circulator);

    ++begin;
    ++begin;
    ++begin;

    EXPECT_TRUE(begin == end);
}

TEST(CirculatorIteratorTest, WorksWithStdAlgorithms) {
    std::vector<int> data = {1, 2, 3, 4, 5};
    MockCirculator<int> circulator(&data, 0);

    auto begin = CirculatorIterator<MockCirculator<int>, int>(circulator);
    auto end = CirculatorIterator<MockCirculator<int>, int>(circulator);

    int sum = 0;
    std::for_each(begin, end, [&sum](int x) { sum += x; });

    EXPECT_EQ(sum, 15);
}

TEST(CirculatorIteratorTest, SingleElementCirculator) {
    std::vector<int> data = {42};
    MockCirculator<int> circulator(&data, 0);

    auto begin = CirculatorIterator<MockCirculator<int>, int>(circulator);
    auto end = CirculatorIterator<MockCirculator<int>, int>(circulator);

    std::vector<int> result;
    for (auto it = begin; it != end; ++it) {
        result.push_back(*it);
    }

    EXPECT_EQ(result.size(), 1);
    EXPECT_EQ(result[0], 42);
}

TEST(CirculatorIteratorTest, StartsFromNonZeroIndex) {
    std::vector<int> data = {1, 2, 3, 4};
    MockCirculator<int> circulator(&data, 2);

    auto begin = CirculatorIterator<MockCirculator<int>, int>(circulator);
    auto end = CirculatorIterator<MockCirculator<int>, int>(circulator);

    std::vector<int> result;
    for (auto it = begin; it != end; ++it) {
        result.push_back(*it);
    }

    std::vector<int> expected = {3, 4, 1, 2};
    EXPECT_EQ(result, expected);
}
