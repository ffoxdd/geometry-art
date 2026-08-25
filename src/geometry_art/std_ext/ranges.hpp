#ifndef GEOMETRY_ART_HELPERS_RANGES_HPP_
#define GEOMETRY_ART_HELPERS_RANGES_HPP_

#include <ranges>
#include <iterator>
#include <algorithm>
#include <concepts>

namespace geometry_art::std_ext {

template<typename Range, typename ValueType>
concept RangeOf =
    std::ranges::range<Range> &&
    std::same_as<std::ranges::range_value_t<Range>, ValueType>;

template<std::ranges::range RangeType>
class CircularPairsView : public std::ranges::view_interface<CircularPairsView<RangeType>> {
 public:
    class Iterator {
     public:
        using iterator_t = std::ranges::iterator_t<const RangeType>;
        using value_type = std::pair<
            std::ranges::range_value_t<const RangeType>,
            std::ranges::range_value_t<const RangeType>
        >;
        using difference_type = std::ptrdiff_t;
        using iterator_category = std::input_iterator_tag;
        using iterator_concept = std::input_iterator_tag;

        Iterator() = default;
        Iterator(iterator_t it, iterator_t begin_it, iterator_t end_it, bool is_wraparound)
            : _it(it), _begin_it(begin_it), _end_it(end_it), _is_wraparound(is_wraparound) {}

        value_type operator*() const {
            if (_is_wraparound) {
                iterator_t last = _begin_it;
                for (auto temp = _begin_it; temp != _end_it; ++temp) {
                    last = temp;
                }
                return {*last, *_begin_it};
            } else {
                auto next = std::next(_it);
                return {*_it, *next};
            }
        }

        Iterator &operator++() {
            if (_is_wraparound) {
                _it = _end_it;
                _is_wraparound = false;
            } else {
                ++_it;
                auto next = std::next(_it);
                if (next == _end_it) {
                    _is_wraparound = true;
                }
            }
            return *this;
        }

        Iterator operator++(int) {
            Iterator tmp = *this;
            ++(*this);
            return tmp;
        }

        bool operator==(const Iterator &other) const {
            return _it == other._it && _is_wraparound == other._is_wraparound;
        }

        bool operator!=(const Iterator &other) const {
            return !(*this == other);
        }

     private:
        iterator_t _it{};
        iterator_t _begin_it{};
        iterator_t _end_it{};
        bool _is_wraparound = false;
    };

    CircularPairsView() = default;
    explicit CircularPairsView(const RangeType &range) : _range_ptr(&range) {}

    Iterator begin() const {
        auto it = std::ranges::begin(*_range_ptr);
        auto end_it = std::ranges::end(*_range_ptr);
        bool is_wraparound = false;

        if (it == end_it) {
            return Iterator(end_it, it, end_it, false);
        }

        auto next = std::next(it);
        if (next == end_it) {
            is_wraparound = true;
        }

        return Iterator(it, it, end_it, is_wraparound);
    }

    Iterator end() const {
        auto end_it = std::ranges::end(*_range_ptr);
        return Iterator(end_it, std::ranges::begin(*_range_ptr), end_it, false);
    }

 private:
    const RangeType *_range_ptr = nullptr;
};

template<std::ranges::range RangeType>
auto circular_adjacent_pairs(const RangeType &range) {
    return CircularPairsView<RangeType>(range);
}

template<std::ranges::range RangeType, typename Predicate>
bool all_circular_adjacent_pairs(const RangeType &range, Predicate pred) {
    auto pairs = circular_adjacent_pairs(range);
    return std::all_of(pairs.begin(), pairs.end(), pred);
}

template<std::ranges::range RangeType>
auto sum(const RangeType &range) {
    std::ranges::range_value_t<RangeType> result{};
    for (const auto &value : range) {
        result += value;
    }
    return result;
}

} // namespace geometry_art::std_ext

namespace geometry_art {
using std_ext::RangeOf;
}

#endif //GEOMETRY_ART_HELPERS_RANGES_HPP_

