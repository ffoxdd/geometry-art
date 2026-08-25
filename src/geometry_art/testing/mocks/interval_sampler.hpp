#ifndef GEOMETRY_ART_TESTING_MOCKS_INTERVAL_SAMPLER_HPP_
#define GEOMETRY_ART_TESTING_MOCKS_INTERVAL_SAMPLER_HPP_

#include "../../math/interval.hpp"
#include <vector>
#include <utility>

namespace geometry_art::testing::mocks {

class MockIntervalSampler {
 public:
    explicit MockIntervalSampler(std::vector<double> sequence)
        : _sequence(std::move(sequence)), _index(0) {}

    [[nodiscard]] double sample(const Interval& interval) {
        double t = _sequence[_index % _sequence.size()];
        _index++;
        return interval.low() + t * (interval.high() - interval.low());
    }

 private:
    std::vector<double> _sequence;
    size_t _index;
};

} // namespace geometry_art::testing::mocks

#endif //GEOMETRY_ART_TESTING_MOCKS_INTERVAL_SAMPLER_HPP_
