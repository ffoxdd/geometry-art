#ifndef GEOMETRY_ART_DLA_SPAWN_SHELL_SPAWN_SHELL_HPP_
#define GEOMETRY_ART_DLA_SPAWN_SHELL_SPAWN_SHELL_HPP_

#include "../../types.hpp"
#include <concepts>

namespace geometry_art::dla {

// The origin-centered sphere through which walkers enter the simulation.
// spawn admits a walker arriving from infinity onto a shell of the given
// radius; confine takes a walker's position after a jump and resolves any
// excursion beyond the shell into the next position to continue from,
// returning positions at or inside the shell untouched.
template<typename T>
concept SpawnShell = requires(T shell, const Vector3& position, double radius) {
    { shell.spawn(radius) } -> std::convertible_to<Vector3>;
    { shell.confine(position, radius) } -> std::convertible_to<Vector3>;
};

} // namespace geometry_art::dla

#endif //GEOMETRY_ART_DLA_SPAWN_SHELL_SPAWN_SHELL_HPP_
