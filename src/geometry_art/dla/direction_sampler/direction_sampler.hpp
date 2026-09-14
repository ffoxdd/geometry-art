#ifndef GEOMETRY_ART_DLA_DIRECTION_SAMPLER_DIRECTION_SAMPLER_HPP_
#define GEOMETRY_ART_DLA_DIRECTION_SAMPLER_DIRECTION_SAMPLER_HPP_

#include "../../types.hpp"
#include <concepts>

namespace geometry_art::dla {

template<typename T>
concept DirectionSampler = requires(T sampler) {
    { sampler.sample() } -> std::convertible_to<Vector3>;
};

} // namespace geometry_art::dla

#endif //GEOMETRY_ART_DLA_DIRECTION_SAMPLER_DIRECTION_SAMPLER_HPP_
