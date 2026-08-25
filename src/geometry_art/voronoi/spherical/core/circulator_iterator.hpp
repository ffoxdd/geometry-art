#ifndef GEOMETRY_ART_VORONOI_SPHERICAL_CORE_CIRCULATOR_ITERATOR_HPP_
#define GEOMETRY_ART_VORONOI_SPHERICAL_CORE_CIRCULATOR_ITERATOR_HPP_

#include <iterator>
#include <cstddef>

namespace geometry_art::voronoi::spherical {

template<typename CirculatorType, typename HandleType>
class CirculatorIterator {
 public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = HandleType;
    using difference_type = std::ptrdiff_t;
    using pointer = const HandleType *;
    using reference = const HandleType &;

    CirculatorIterator() :
        _circulator(),
        _start(true) {
    }

    explicit CirculatorIterator(CirculatorType circulator) :
        _circulator(circulator),
        _start(true) {
    }

    reference operator*() const {
        return *_circulator;;
    }

    pointer operator->() const {
        return &(*_circulator);
    }

    CirculatorIterator &operator++() {
        ++_circulator;
        _start = false;
        return *this;
    }

    CirculatorIterator operator++(int) {
        CirculatorIterator tmp = *this;
        ++(*this);
        return tmp;
    }

    friend bool operator==(const CirculatorIterator &a, const CirculatorIterator &b) {
        if (a._circulator == b._circulator && a._start && b._start) {
            return false;
        }

        return a._circulator == b._circulator;
    }

    friend bool operator!=(const CirculatorIterator &a, const CirculatorIterator &b) {
        return !(a == b);
    }

 private:
    CirculatorType _circulator;
    bool _start;
};

} // namespace geometry_art::voronoi::spherical

#endif //GEOMETRY_ART_VORONOI_SPHERICAL_CORE_CIRCULATOR_ITERATOR_HPP_

