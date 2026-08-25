#ifndef GEOMETRY_ART_MATH_CIRCULAR_INTERVAL_SAMPLER_UNIFORM_CIRCULAR_INTERVAL_SAMPLER_HPP_
#define GEOMETRY_ART_MATH_CIRCULAR_INTERVAL_SAMPLER_UNIFORM_CIRCULAR_INTERVAL_SAMPLER_HPP_

#include "../circular_interval.hpp"
#include <random>

namespace geometry_art::math {

class UniformCircularIntervalSampler {
 public:
    UniformCircularIntervalSampler() : _random_engine(std::random_device{}()) {}

    explicit UniformCircularIntervalSampler(unsigned int seed) : _random_engine(seed) {}

    template<double PERIOD>
    [[nodiscard]] inline double sample(const CircularInterval<PERIOD>& interval) {
        double u = _distribution(_random_engine);
        return interval.start() + u * interval.measure();
    }

 private:
    std::mt19937 _random_engine;
    std::uniform_real_distribution<double> _distribution{0.0, 1.0};
};

} // namespace geometry_art::math

namespace geometry_art {
using UniformCircularIntervalSampler = math::UniformCircularIntervalSampler;
}

#endif //GEOMETRY_ART_MATH_CIRCULAR_INTERVAL_SAMPLER_UNIFORM_CIRCULAR_INTERVAL_SAMPLER_HPP_
