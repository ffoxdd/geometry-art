#ifndef GEOMETRY_ART_TESTING_MOCKS_DIRECTION_SAMPLER_HPP_
#define GEOMETRY_ART_TESTING_MOCKS_DIRECTION_SAMPLER_HPP_

#include "../../types.hpp"
#include <vector>
#include <utility>

namespace geometry_art::testing::mocks {

class MockDirectionSampler {
 public:
    explicit MockDirectionSampler(std::vector<Vector3> sequence)
        : _sequence(std::move(sequence)), _index(0) {}

    [[nodiscard]] Vector3 sample() {
        Vector3 direction = _sequence[_index % _sequence.size()];
        _index++;
        return direction;
    }

 private:
    std::vector<Vector3> _sequence;
    size_t _index;
};

} // namespace geometry_art::testing::mocks

#endif //GEOMETRY_ART_TESTING_MOCKS_DIRECTION_SAMPLER_HPP_
