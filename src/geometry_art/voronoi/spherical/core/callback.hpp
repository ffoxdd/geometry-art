#ifndef GEOMETRY_ART_VORONOI_SPHERICAL_CORE_CALLBACK_HPP_
#define GEOMETRY_ART_VORONOI_SPHERICAL_CORE_CALLBACK_HPP_

#include "sphere.hpp"
#include <functional>

namespace geometry_art::voronoi::spherical {

using Callback = std::function<void(const Sphere&)>;

inline Callback noop_callback() {
    return [](const Sphere&) {};
}

} // namespace geometry_art::voronoi::spherical

#endif //GEOMETRY_ART_VORONOI_SPHERICAL_CORE_CALLBACK_HPP_
