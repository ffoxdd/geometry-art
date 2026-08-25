#ifndef GEOMETRY_ART_MATH_INTERVAL_SAMPLER_INTERVAL_SAMPLER_HPP_
#define GEOMETRY_ART_MATH_INTERVAL_SAMPLER_INTERVAL_SAMPLER_HPP_

#include "../interval.hpp"
#include <concepts>

namespace geometry_art::math {

template<typename T>
concept IntervalSampler = requires(T sampler, const Interval& interval) {
    { sampler.sample(interval) } -> std::convertible_to<double>;
};

} // namespace geometry_art::math

namespace geometry_art {
using math::IntervalSampler;
}

#endif //GEOMETRY_ART_MATH_INTERVAL_SAMPLER_INTERVAL_SAMPLER_HPP_
