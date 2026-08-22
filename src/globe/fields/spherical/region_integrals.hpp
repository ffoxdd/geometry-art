#ifndef GLOBEART_SRC_GLOBE_FIELDS_SPHERICAL_REGION_INTEGRALS_HPP_
#define GLOBEART_SRC_GLOBE_FIELDS_SPHERICAL_REGION_INTEGRALS_HPP_

#include "../../types.hpp"

namespace globe::fields::spherical {

struct RegionIntegrals {
    double mass;
    Vector3 first_moment;
};

} // namespace globe::fields::spherical

#endif //GLOBEART_SRC_GLOBE_FIELDS_SPHERICAL_REGION_INTEGRALS_HPP_
