#ifndef GEOMETRY_ART_GEOMETRY_CARTESIAN_BOUNDING_BOX_HPP_
#define GEOMETRY_ART_GEOMETRY_CARTESIAN_BOUNDING_BOX_HPP_

#include "../../math/interval.hpp"

namespace geometry_art::geometry::cartesian {

using geometry_art::Interval;

class BoundingBox {
 public:
    BoundingBox();
    BoundingBox(Interval x_interval, Interval y_interval, Interval z_interval);

    [[nodiscard]] Interval x_interval() const;
    [[nodiscard]] Interval y_interval() const;
    [[nodiscard]] Interval z_interval() const;

    [[nodiscard]] static BoundingBox unit_cube();

 private:
    const Interval _x_interval;
    const Interval _y_interval;
    const Interval _z_interval;
};

inline BoundingBox::BoundingBox() :
    BoundingBox(
        Interval(-1, 1),
        Interval(-1, 1),
        Interval(-1, 1)
    ) {
}

inline BoundingBox::BoundingBox(Interval x_interval, Interval y_interval, Interval z_interval) :
    _x_interval(x_interval),
    _y_interval(y_interval),
    _z_interval(z_interval) {
}

inline Interval BoundingBox::x_interval() const {
    return _x_interval;
}

inline Interval BoundingBox::y_interval() const {
    return _y_interval;
}

inline Interval BoundingBox::z_interval() const {
    return _z_interval;
}

inline BoundingBox BoundingBox::unit_cube() {
    return BoundingBox(Interval(-1, 1), Interval(-1, 1), Interval(-1, 1));
}

} // namespace geometry_art::geometry::cartesian

namespace geometry_art {
using geometry::cartesian::BoundingBox;
}

#endif //GEOMETRY_ART_GEOMETRY_CARTESIAN_BOUNDING_BOX_HPP_
