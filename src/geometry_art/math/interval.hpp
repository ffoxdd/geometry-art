#ifndef GEOMETRY_ART_NOISE_GENERATOR_RANGE_H_
#define GEOMETRY_ART_NOISE_GENERATOR_RANGE_H_

#include "../std_ext/ranges.hpp"
#include <CGAL/assertions.h>
#include <limits>
#include <concepts>
#include <ranges>
#include <algorithm>
#include <numeric>

namespace geometry_art::math {

class Interval {
 public:
    Interval();
    constexpr Interval(double low, double high);
    template<RangeOf<double> DoubleRange> explicit Interval(const DoubleRange &range);

    [[nodiscard]] constexpr double low() const;
    [[nodiscard]] constexpr double high() const;
    [[nodiscard]] constexpr double measure() const;
    [[nodiscard]] constexpr double midpoint() const;
    [[nodiscard]] constexpr double clamp(double value) const;
    [[nodiscard]] constexpr double normalize(double value) const;
    [[nodiscard]] constexpr bool contains(double value) const;

    [[nodiscard]] static constexpr Interval hull(const Interval &a, const Interval &b);
    [[nodiscard]] static constexpr double remap(double value, const Interval& source, const Interval& target);
    [[nodiscard]] static constexpr Interval hull(const Interval &interval, double value);
    [[nodiscard]] static constexpr Interval hull(double a, double b);

    template<RangeOf<Interval> IntervalRange>
    [[nodiscard]] static Interval hull(const IntervalRange &intervals);

 private:
    double _low;
    double _high;
};

inline Interval::Interval() :
    _low(std::numeric_limits<double>::max()),
    _high(std::numeric_limits<double>::lowest()) {
};

constexpr Interval::Interval(double low, double high) :
    _low(low),
    _high(high) {
}

inline const Interval UNIT_INTERVAL{0.0, 1.0};

template<RangeOf<double> DoubleRange>
Interval::Interval(const DoubleRange &range) {
    auto [min_it, max_it] = std::ranges::minmax_element(range);

    _low = *min_it;
    _high = *max_it;
}

constexpr double Interval::low() const {
    return _low;
}

constexpr double Interval::high() const {
    return _high;
}

constexpr double Interval::measure() const {
    return _high - _low;
}

constexpr double Interval::midpoint() const {
    return (_low + _high) / 2.0;
}

constexpr double Interval::clamp(double value) const {
    return std::clamp(value, _low, _high);
}

constexpr double Interval::normalize(double value) const {
    return remap(value, *this, UNIT_INTERVAL);
}

constexpr double Interval::remap(double value, const Interval& source, const Interval& target) {
    double t = (value - source.low()) / source.measure();
    t = std::clamp(t, 0.0, 1.0);
    return target.low() + t * target.measure();
}

constexpr bool Interval::contains(double value) const {
    return _low <= value && value <= _high;
}

constexpr Interval Interval::hull(const Interval &a, const Interval &b) {
    return Interval(
        std::min(a._low, b._low),
        std::max(a._high, b._high)
    );
}

constexpr Interval Interval::hull(const Interval &interval, double value) {
    return Interval(
        std::min(interval._low, value),
        std::max(interval._high, value)
    );
}

constexpr Interval Interval::hull(double a, double b) {
    return Interval(
        std::min(a, b),
        std::max(a, b)
    );
}

template<RangeOf<Interval> IntervalRange>
Interval Interval::hull(const IntervalRange &intervals) {
    CGAL_precondition(!std::ranges::empty(intervals));

    return std::accumulate(
        std::next(std::ranges::begin(intervals)),
        std::ranges::end(intervals),
        *std::ranges::begin(intervals),
        [](const Interval &a, const Interval &b) { return hull(a, b); }
    );
}

} // namespace geometry_art::math

namespace geometry_art {
using Interval = math::Interval;
inline const Interval UNIT_INTERVAL = math::UNIT_INTERVAL;
}

#endif //GEOMETRY_ART_NOISE_GENERATOR_RANGE_H_
