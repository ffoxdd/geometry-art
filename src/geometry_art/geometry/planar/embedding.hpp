#ifndef GEOMETRY_ART_GEOMETRY_PLANAR_EMBEDDING_HPP_
#define GEOMETRY_ART_GEOMETRY_PLANAR_EMBEDDING_HPP_

#include "domain.hpp"
#include "../../types.hpp"
#include <algorithm>
#include <cmath>

namespace geometry_art::geometry::planar {

// Where a flat domain lives in space, read off how it closes: the plane as
// itself, the cylinder rolled about the vertical axis, the torus as a ring
// whose tube is the wrapped height. Each is oriented so that a
// counter-clockwise turn in the rectangle is counter-clockwise seen from
// outside, which is the side the normal points to.
class Embedding {
 public:
    [[nodiscard]] static Embedding of(const Domain& domain);

    [[nodiscard]] Vector3 position(const Vector2& point) const;
    [[nodiscard]] Vector3 normal(const Vector2& point) const;
    [[nodiscard]] bool curved() const { return _domain.has_wrapped_axis(); }

 private:
    static constexpr double LARGEST_TUBE_SHARE = 0.75;

    Domain _domain;

    explicit Embedding(const Domain& domain);

    [[nodiscard]] double radius() const { return _domain.width / TWO_PI; }
    [[nodiscard]] double tube() const { return std::min(_domain.height / TWO_PI, LARGEST_TUBE_SHARE * radius()); }
    [[nodiscard]] double around(const Vector2& point) const { return TWO_PI * point.x() / _domain.width; }
    [[nodiscard]] double along(const Vector2& point) const { return TWO_PI * point.y() / _domain.height; }
};

inline Embedding Embedding::of(const Domain& domain) {
    return Embedding(domain);
}

inline Embedding::Embedding(const Domain& domain) :
    _domain(domain) {
}

inline Vector3 Embedding::position(const Vector2& point) const {
    if (!_domain.wrapped(0)) {
        return Vector3(point.x() - 0.5 * _domain.width, point.y() - 0.5 * _domain.height, 0.0);
    }

    double turn = around(point);

    if (!_domain.wrapped(1)) {
        return Vector3(radius() * std::cos(turn), point.y() - 0.5 * _domain.height, -radius() * std::sin(turn));
    }

    double lift = along(point);
    double ring = radius() + tube() * std::cos(lift);

    return Vector3(ring * std::cos(turn), tube() * std::sin(lift), -ring * std::sin(turn));
}

inline Vector3 Embedding::normal(const Vector2& point) const {
    if (!_domain.wrapped(0)) {
        return Vector3(0.0, 0.0, 1.0);
    }

    double turn = around(point);

    if (!_domain.wrapped(1)) {
        return Vector3(std::cos(turn), 0.0, -std::sin(turn));
    }

    double lift = along(point);

    return Vector3(std::cos(lift) * std::cos(turn), std::sin(lift), -std::cos(lift) * std::sin(turn));
}

} // namespace geometry_art::geometry::planar

#endif //GEOMETRY_ART_GEOMETRY_PLANAR_EMBEDDING_HPP_
