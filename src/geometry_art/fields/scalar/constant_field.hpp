#ifndef GEOMETRY_ART_FIELDS_SCALAR_CONSTANT_FIELD_HPP_
#define GEOMETRY_ART_FIELDS_SCALAR_CONSTANT_FIELD_HPP_

#include "../../types.hpp"

namespace geometry_art::fields::scalar {

class ConstantField {
public:
    explicit ConstantField(double value = 1.0) :
        _value(value) {}

    inline double value(const VectorS2 &) const {
        return _value;
    }

    inline double max_frequency() const {
        return 0.0;
    }

private:
    double _value;
};

} // namespace geometry_art::fields::scalar

#endif //GEOMETRY_ART_FIELDS_SCALAR_CONSTANT_FIELD_HPP_
