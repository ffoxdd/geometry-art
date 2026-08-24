#ifndef GLOBEART_SRC_GLOBE_FIELDS_REGION_INTEGRALS_HPP_
#define GLOBEART_SRC_GLOBE_FIELDS_REGION_INTEGRALS_HPP_

#include "../types.hpp"

namespace globe::fields {

struct RegionIntegrals {
    double mass;
    Vector3 first_moment;
};

} // namespace globe::fields

#endif //GLOBEART_SRC_GLOBE_FIELDS_REGION_INTEGRALS_HPP_
