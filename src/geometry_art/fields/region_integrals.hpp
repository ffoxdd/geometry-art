#ifndef GEOMETRY_ART_FIELDS_REGION_INTEGRALS_HPP_
#define GEOMETRY_ART_FIELDS_REGION_INTEGRALS_HPP_

#include "../types.hpp"

namespace geometry_art::fields {

struct RegionIntegrals {
    double mass;
    Vector3 first_moment;
};

} // namespace geometry_art::fields

#endif //GEOMETRY_ART_FIELDS_REGION_INTEGRALS_HPP_
