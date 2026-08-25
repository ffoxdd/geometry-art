#ifndef GEOMETRY_ART_MATH_CIRCULAR_INTERVAL_HPP_
#define GEOMETRY_ART_MATH_CIRCULAR_INTERVAL_HPP_

#include "../std_ext/ranges.hpp"
#include <CGAL/assertions.h>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <ranges>

namespace geometry_art::math {

template<double PERIOD>
class CircularInterval {
 public:
    CircularInterval(double start, double measure);

    static CircularInterval from_to(double from, double to);
    static CircularInterval full();

    [[nodiscard]] double start() const;
    [[nodiscard]] double measure() const;
    [[nodiscard]] double end() const;
    [[nodiscard]] bool is_full() const;
    [[nodiscard]] bool contains(double value) const;

    [[nodiscard]] static CircularInterval hull(
        const CircularInterval &a,
        const CircularInterval &b
    );

    template<RangeOf<CircularInterval> CircularIntervalRange>
    [[nodiscard]] static CircularInterval hull(const CircularIntervalRange &intervals);

 private:
    double _start;
    double _measure;

    [[nodiscard]] static double normalize(double value);
    [[nodiscard]] static double forward_distance(double from, double to);

    [[nodiscard]] static double measure_to_cover(
        double start,
        double base_measure,
        const CircularInterval &other
    );
};

template<double PERIOD>
CircularInterval<PERIOD>::CircularInterval(double start, double measure) :
    _start(normalize(start)),
    _measure(std::clamp(measure, 0.0, PERIOD)) {
}

template<double PERIOD>
CircularInterval<PERIOD> CircularInterval<PERIOD>::from_to(double from, double to) {
    double start = normalize(from);
    double measure = forward_distance(start, normalize(to));

    if (measure == 0.0) {
        measure = PERIOD;
    }

    return CircularInterval(start, measure);
}

template<double PERIOD>
CircularInterval<PERIOD> CircularInterval<PERIOD>::full() {
    return CircularInterval(0.0, PERIOD);
}

template<double PERIOD>
double CircularInterval<PERIOD>::start() const {
    return _start;
}

template<double PERIOD>
double CircularInterval<PERIOD>::measure() const {
    return _measure;
}

template<double PERIOD>
double CircularInterval<PERIOD>::end() const {
    return _start + _measure;
}

template<double PERIOD>
bool CircularInterval<PERIOD>::is_full() const {
    return _measure >= PERIOD;
}

template<double PERIOD>
bool CircularInterval<PERIOD>::contains(double value) const {
    if (is_full()) {
        return true;
    }

    return forward_distance(_start, normalize(value)) <= _measure;
}

template<double PERIOD>
CircularInterval<PERIOD> CircularInterval<PERIOD>::hull(
    const CircularInterval &a,
    const CircularInterval &b
) {
    if (a.is_full() || b.is_full()) {
        return full();
    }

    double measure_from_a = measure_to_cover(a._start, a._measure, b);
    double measure_from_b = measure_to_cover(b._start, b._measure, a);

    if (measure_from_a >= PERIOD && measure_from_b >= PERIOD) {
        return full();
    }

    if (measure_from_a <= measure_from_b) {
        return CircularInterval(a._start, measure_from_a);
    } else {
        return CircularInterval(b._start, measure_from_b);
    }
}

template<double PERIOD>
template<RangeOf<CircularInterval<PERIOD>> CircularIntervalRange>
CircularInterval<PERIOD> CircularInterval<PERIOD>::hull(const CircularIntervalRange &intervals) {
    CGAL_precondition(!std::ranges::empty(intervals));

    return std::accumulate(
        std::next(std::ranges::begin(intervals)),
        std::ranges::end(intervals),
        *std::ranges::begin(intervals),
        [](const CircularInterval &a, const CircularInterval &b) { return hull(a, b); }
    );
}

template<double PERIOD>
double CircularInterval<PERIOD>::normalize(double value) {
    double result = std::fmod(value, PERIOD);

    if (result < 0.0) {
        result += PERIOD;
    }

    return result;
}

template<double PERIOD>
double CircularInterval<PERIOD>::forward_distance(double from, double to) {
    double distance = to - from;

    if (distance < 0.0) {
        distance += PERIOD;
    }

    return distance;
}

template<double PERIOD>
double CircularInterval<PERIOD>::measure_to_cover(
    double start,
    double base_measure,
    const CircularInterval &other
) {
    if (other._measure >= PERIOD) {
        return PERIOD;
    }

    double dist_to_other_start = forward_distance(start, other._start);
    double dist_to_other_end = dist_to_other_start + other._measure;

    if (dist_to_other_end <= base_measure) {
        return base_measure;
    }

    if (dist_to_other_start <= base_measure) {
        return std::min(dist_to_other_end, PERIOD);
    }

    if (dist_to_other_end > PERIOD) {
        double wrapped_end = dist_to_other_end - PERIOD;

        if (wrapped_end <= base_measure) {
            return PERIOD;
        }
    }

    return std::min(dist_to_other_end, PERIOD);
}

} // namespace geometry_art::math

namespace geometry_art {
template<double PERIOD>
using CircularInterval = math::CircularInterval<PERIOD>;
}

#endif //GEOMETRY_ART_MATH_CIRCULAR_INTERVAL_HPP_
