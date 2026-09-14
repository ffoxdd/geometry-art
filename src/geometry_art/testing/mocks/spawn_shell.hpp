#ifndef GEOMETRY_ART_TESTING_MOCKS_SPAWN_SHELL_HPP_
#define GEOMETRY_ART_TESTING_MOCKS_SPAWN_SHELL_HPP_

#include "../../types.hpp"
#include <utility>

namespace geometry_art::testing::mocks {

// Spawns at a scripted absolute position and never intervenes in the walk.
class MockSpawnShell {
 public:
    explicit MockSpawnShell(Vector3 spawn_position)
        : _spawn_position(std::move(spawn_position)) {}

    [[nodiscard]] Vector3 spawn(double) const {
        return _spawn_position;
    }

    [[nodiscard]] Vector3 confine(const Vector3& position, double) const {
        return position;
    }

 private:
    Vector3 _spawn_position;
};

} // namespace geometry_art::testing::mocks

#endif //GEOMETRY_ART_TESTING_MOCKS_SPAWN_SHELL_HPP_
