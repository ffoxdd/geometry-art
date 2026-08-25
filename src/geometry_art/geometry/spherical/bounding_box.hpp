#ifndef GEOMETRY_ART_GEOMETRY_SPHERICAL_BOUNDING_BOX_HPP_
#define GEOMETRY_ART_GEOMETRY_SPHERICAL_BOUNDING_BOX_HPP_

#include "../../math/interval.hpp"
#include "../../math/circular_interval.hpp"
#include "../../types.hpp"
#include <cmath>

namespace geometry_art::geometry::spherical {

using geometry_art::Interval;
using geometry_art::CircularInterval;
using geometry_art::TWO_PI;
using geometry_art::VectorS2;

using ThetaInterval = CircularInterval<TWO_PI>;

class BoundingBox {
 public:
    BoundingBox(ThetaInterval theta_interval, Interval z_interval);

    [[nodiscard]] ThetaInterval theta_interval() const;
    [[nodiscard]] Interval z_interval() const;
    [[nodiscard]] double area() const;
    [[nodiscard]] VectorS2 center() const;
    [[nodiscard]] double bounding_sphere_radius() const;
    [[nodiscard]] bool contains(const VectorS2& point) const;

    [[nodiscard]] static BoundingBox full_sphere();

 private:
    ThetaInterval _theta_interval;
    Interval _z_interval;

    [[nodiscard]] static double theta(double x, double y);
};

inline BoundingBox::BoundingBox(ThetaInterval theta_interval, Interval z_interval) :
    _theta_interval(theta_interval),
    _z_interval(z_interval) {
}

inline ThetaInterval BoundingBox::theta_interval() const {
    return _theta_interval;
}

inline Interval BoundingBox::z_interval() const {
    return _z_interval;
}

inline double BoundingBox::area() const {
    return _theta_interval.measure() * _z_interval.measure();
}

inline VectorS2 BoundingBox::center() const {
    double theta_mid = _theta_interval.start() + _theta_interval.measure() / 2.0;
    if (theta_mid >= TWO_PI) {
        theta_mid -= TWO_PI;
    }

    double z_mid = _z_interval.midpoint();
    double r_mid = std::sqrt(1.0 - z_mid * z_mid);

    return VectorS2(
        r_mid * std::cos(theta_mid),
        r_mid * std::sin(theta_mid),
        z_mid
    );
}

inline double BoundingBox::bounding_sphere_radius() const {
    double z_span = _z_interval.measure();
    double theta_span = _theta_interval.measure();

    double r_low = std::sqrt(std::max(0.0, 1.0 - _z_interval.low() * _z_interval.low()));
    double r_high = std::sqrt(std::max(0.0, 1.0 - _z_interval.high() * _z_interval.high()));
    double r_max = std::max(r_low, r_high);
    double chord = 2.0 * r_max * std::sin(theta_span / 2.0);

    return std::sqrt(z_span * z_span + chord * chord);
}

inline bool BoundingBox::contains(const VectorS2& point) const {
    double theta_val = theta(point.x(), point.y());
    double z_val = point.z();
    return _theta_interval.contains(theta_val) && _z_interval.contains(z_val);
}

inline BoundingBox BoundingBox::full_sphere() {
    return BoundingBox(ThetaInterval::full(), Interval(-1, 1));
}

inline double BoundingBox::theta(double x, double y) {
    double t = std::atan2(y, x);
    return t < 0.0 ? t + TWO_PI : t;
}

} // namespace geometry_art::geometry::spherical

namespace geometry_art {
using SphericalBoundingBox = geometry::spherical::BoundingBox;
using ThetaInterval = geometry::spherical::ThetaInterval;
}

#endif //GEOMETRY_ART_GEOMETRY_SPHERICAL_BOUNDING_BOX_HPP_
