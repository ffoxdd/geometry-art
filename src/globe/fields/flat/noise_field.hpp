#ifndef GLOBEART_SRC_GLOBE_FIELDS_FLAT_NOISE_FIELD_HPP_
#define GLOBEART_SRC_GLOBE_FIELDS_FLAT_NOISE_FIELD_HPP_

#include "../scalar/noise_field.hpp"
#include "../../math/interval.hpp"
#include "../../types.hpp"
#include <cmath>

namespace globe::fields::flat {

using globe::math::Interval;

// Periodic noise on a rectangle: the rectangle is embedded as a torus in
// space and the same solid noise the sphere samples is read along it, so
// wrapping in either direction is seamless by construction.
class NoiseField {
 public:
    NoiseField(double width, double height, Interval output_range = Interval(0, 1), int seed = 1546);

    [[nodiscard]] double value(const Vector2& point);

 private:
    static constexpr double MAJOR_RADIUS = 1.0;
    static constexpr double MINOR_RADIUS = 0.45;

    double _width;
    double _height;
    scalar::NoiseField _noise;
};

inline NoiseField::NoiseField(double width, double height, Interval output_range, int seed) :
    _width(width),
    _height(height),
    _noise(output_range, seed) {
}

inline double NoiseField::value(const Vector2& point) {
    double around = 2.0 * M_PI * point.x() / _width;
    double along = 2.0 * M_PI * point.y() / _height;
    double ring = MAJOR_RADIUS + MINOR_RADIUS * std::cos(along);

    return _noise.value(VectorS2(
        ring * std::cos(around),
        ring * std::sin(around),
        MINOR_RADIUS * std::sin(along)
    ));
}

} // namespace globe::fields::flat

#endif //GLOBEART_SRC_GLOBE_FIELDS_FLAT_NOISE_FIELD_HPP_
