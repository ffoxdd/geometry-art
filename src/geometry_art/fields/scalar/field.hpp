#ifndef GEOMETRY_ART_FIELDS_SCALAR_FIELD_HPP_
#define GEOMETRY_ART_FIELDS_SCALAR_FIELD_HPP_

#include "../../types.hpp"

namespace geometry_art::fields::scalar {

template<typename T>
concept Field = requires( T field, const VectorS2 &point ) {
    { field.value(point) } -> std::convertible_to<double>;
};

} // namespace geometry_art::fields::scalar

#endif //GEOMETRY_ART_FIELDS_SCALAR_FIELD_HPP_
