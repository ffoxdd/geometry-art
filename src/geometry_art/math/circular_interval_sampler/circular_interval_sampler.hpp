#ifndef GEOMETRY_ART_MATH_CIRCULAR_INTERVAL_SAMPLER_CIRCULAR_INTERVAL_SAMPLER_HPP_
#define GEOMETRY_ART_MATH_CIRCULAR_INTERVAL_SAMPLER_CIRCULAR_INTERVAL_SAMPLER_HPP_

#include "../circular_interval.hpp"
#include <concepts>

namespace geometry_art::math {

template<typename T, double PERIOD>
concept CircularIntervalSampler = requires(T sampler, const CircularInterval<PERIOD>& interval) {
    { sampler.sample(interval) } -> std::convertible_to<double>;
};

} // namespace geometry_art::math

namespace geometry_art {
using math::CircularIntervalSampler;
}

#endif //GEOMETRY_ART_MATH_CIRCULAR_INTERVAL_SAMPLER_CIRCULAR_INTERVAL_SAMPLER_HPP_
