#ifndef GEOMETRY_ART_MATH_INTERVAL_SAMPLER_UNIFORM_INTERVAL_SAMPLER_HPP_
#define GEOMETRY_ART_MATH_INTERVAL_SAMPLER_UNIFORM_INTERVAL_SAMPLER_HPP_

#include "../interval.hpp"
#include <random>

namespace geometry_art::math {

class UniformIntervalSampler {
 public:
    UniformIntervalSampler() : _random_engine(std::random_device{}()) {}

    explicit UniformIntervalSampler(unsigned int seed) : _random_engine(seed) {}

    [[nodiscard]] inline double sample(const Interval& interval) {
        return distribution(interval)(_random_engine);
    }

 private:
    std::mt19937 _random_engine;

    static inline std::uniform_real_distribution<double> distribution(const Interval& interval) {
        return std::uniform_real_distribution<>(interval.low(), interval.high());
    }
};

} // namespace geometry_art::math

namespace geometry_art {
using UniformIntervalSampler = math::UniformIntervalSampler;
}

#endif //GEOMETRY_ART_MATH_INTERVAL_SAMPLER_UNIFORM_INTERVAL_SAMPLER_HPP_
