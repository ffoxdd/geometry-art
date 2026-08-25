#ifndef GEOMETRY_ART_TESTING_MOCKS_CIRCULAR_INTERVAL_SAMPLER_HPP_
#define GEOMETRY_ART_TESTING_MOCKS_CIRCULAR_INTERVAL_SAMPLER_HPP_

#include "../../math/circular_interval.hpp"
#include <vector>
#include <utility>

namespace geometry_art::testing::mocks {

class MockCircularIntervalSampler {
 public:
    explicit MockCircularIntervalSampler(std::vector<double> sequence)
        : _sequence(std::move(sequence)), _index(0) {}

    template<double PERIOD>
    [[nodiscard]] double sample(const CircularInterval<PERIOD>& interval) {
        double t = _sequence[_index % _sequence.size()];
        _index++;
        return interval.start() + t * interval.measure();
    }

 private:
    std::vector<double> _sequence;
    size_t _index;
};

} // namespace geometry_art::testing::mocks

#endif //GEOMETRY_ART_TESTING_MOCKS_CIRCULAR_INTERVAL_SAMPLER_HPP_
