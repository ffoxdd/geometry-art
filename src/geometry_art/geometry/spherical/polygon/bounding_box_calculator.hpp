#ifndef GEOMETRY_ART_GEOMETRY_SPHERICAL_POLYGON_BOUNDING_BOX_CALCULATOR_HPP_
#define GEOMETRY_ART_GEOMETRY_SPHERICAL_POLYGON_BOUNDING_BOX_CALCULATOR_HPP_

// IWYU pragma: private; include "polygon.hpp"

#include "../bounding_box.hpp"
#include "../arc.hpp"
#include "../helpers.hpp"
#include "../../../types.hpp"
#include "../../../math/interval.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <ranges>

namespace geometry_art::geometry::spherical::polygon {

using geometry_art::Interval;
using geometry_art::VectorS2;
using geometry_art::geometry::spherical::Arc;
using geometry_art::geometry::spherical::BoundingBox;
using geometry_art::geometry::spherical::ThetaInterval;

class BoundingBoxCalculator {
 public:
    explicit BoundingBoxCalculator(const std::vector<Arc>& arcs);

    [[nodiscard]] BoundingBox calculate() const;

 private:

    const std::vector<Arc>& _arcs;

    [[nodiscard]] ThetaInterval theta_interval() const;
    [[nodiscard]] Interval z_interval() const;
    [[nodiscard]] auto arc_z_intervals() const;
    [[nodiscard]] auto arc_theta_intervals() const;
    [[nodiscard]] bool contains_point(const VectorS2& point) const;

    [[nodiscard]] static ThetaInterval arc_theta_interval(const Arc& arc);
    [[nodiscard]] static Interval arc_z_interval(const Arc& arc);
    [[nodiscard]] static Interval z_interval_from_endpoints(const Arc& arc);
    [[nodiscard]] static double theta(const VectorS2& point);
    [[nodiscard]] static bool is_zero(const VectorS2& v);
    [[nodiscard]] static bool is_zero(double v);
    [[nodiscard]] static VectorS2 project_to_xy_plane(const VectorS2& v);
    [[nodiscard]] static VectorS2 perpendicular_in_xy_plane(const VectorS2& v);
    [[nodiscard]] static bool parallel(const VectorS2& v1, const VectorS2& v2);
};

inline BoundingBoxCalculator::BoundingBoxCalculator(
    const std::vector<Arc>& arcs
) : _arcs(arcs) {}

inline auto BoundingBoxCalculator::arc_z_intervals() const {
    return _arcs | std::views::transform(
        [](const Arc& arc) { return arc_z_interval(arc); }
    );
}

inline auto BoundingBoxCalculator::arc_theta_intervals() const {
    return _arcs | std::views::transform(
        [](const Arc& arc) { return arc_theta_interval(arc); }
    );
}

inline ThetaInterval BoundingBoxCalculator::theta_interval() const {
    return ThetaInterval::hull(arc_theta_intervals());
}

inline Interval BoundingBoxCalculator::z_interval() const {
    Interval interval = Interval::hull(arc_z_intervals());

    VectorS2 north_pole(0, 0, 1);
    VectorS2 south_pole(0, 0, -1);

    if (contains_point(north_pole)) {
        interval = Interval::hull(interval, north_pole.z());
    }

    if (contains_point(south_pole)) {
        interval = Interval::hull(interval, south_pole.z());
    }

    return interval;
}

inline BoundingBox BoundingBoxCalculator::calculate() const {
    if (_arcs.empty()) {
        return BoundingBox::full_sphere();
    }

    return {theta_interval(), z_interval()};
}

inline bool BoundingBoxCalculator::contains_point(const VectorS2& point) const {
    double angle_sum = 0;

    for (const auto& arc : _arcs) {
        const VectorS2& a = arc.source();
        const VectorS2& b = arc.target();

        VectorS2 ap = a.cross(point).normalized();
        VectorS2 bp = b.cross(point).normalized();

        double angle = distance(ap, bp);
        const VectorS2& normal = arc.normal();

        double sign = normal.dot(point);

        if (is_zero(sign)) {
            return arc.contains(point);
        }

        angle_sum += (sign > 0) ? angle : -angle;
    }

    return std::fabs(angle_sum - (2 * M_PI)) < M_PI;
}

inline Interval BoundingBoxCalculator::arc_z_interval(const Arc& arc) {
    const VectorS2& n = arc.normal();
    Interval z_int = z_interval_from_endpoints(arc);
    VectorS2 n_cross_z = n.cross(VectorS2(0, 0, 1));

    if (is_zero(n_cross_z)) {
        return z_int;
    }

    VectorS2 critical_dir = n_cross_z.cross(n);

    if (is_zero(critical_dir)) {
        return z_int;
    }

    VectorS2 critical_point_up = critical_dir.normalized();
    VectorS2 critical_point_down = -critical_point_up;

    if (arc.contains(critical_point_up)) {
        z_int = Interval::hull(z_int, critical_point_up.z());
    }

    if (arc.contains(critical_point_down)) {
        z_int = Interval::hull(z_int, critical_point_down.z());
    }

    return z_int;
}

inline ThetaInterval BoundingBoxCalculator::arc_theta_interval(const Arc& arc) {
    const VectorS2& n = arc.normal();
    const VectorS2& source_v = arc.source();
    const VectorS2& target_v = arc.target();
    double theta_source = theta(source_v);
    double theta_target = theta(target_v);
    bool is_horizontal = parallel(n, VectorS2(0, 0, 1));

    if (is_horizontal) {
        bool ccw = n.z() > 0;
        if (ccw) {
            return ThetaInterval::from_to(theta_source, theta_target);
        } else {
            return ThetaInterval::from_to(theta_target, theta_source);
        }
    }

    VectorS2 n_xy = project_to_xy_plane(n);

    if (is_zero(n_xy)) {
        return ThetaInterval::full();
    }

    VectorS2 midpoint_dir = source_v + target_v;
    if (is_zero(midpoint_dir)) {
        return ThetaInterval::full();
    }
    VectorS2 midpoint = midpoint_dir.normalized();

    bool short_arc = arc.contains(midpoint);

    double min_theta = std::min(theta_source, theta_target);
    double max_theta = std::max(theta_source, theta_target);
    bool gap_small = (max_theta - min_theta) <= M_PI;

    ThetaInterval result = [&]() {
        if (short_arc) {
            if (gap_small) {
                return ThetaInterval(min_theta, max_theta - min_theta);
            } else {
                return ThetaInterval::from_to(max_theta, min_theta);
            }
        } else {
            if (gap_small) {
                return ThetaInterval::from_to(max_theta, min_theta);
            } else {
                return ThetaInterval(min_theta, max_theta - min_theta);
            }
        }
    }();

    VectorS2 n_xy_normalized = n_xy.normalized();
    VectorS2 radial_dir1 = perpendicular_in_xy_plane(n_xy_normalized);
    VectorS2 radial_dir2 = -radial_dir1;

    VectorS2 critical_dir1 = n.cross(radial_dir1);
    VectorS2 critical_dir2 = n.cross(radial_dir2);

    if (is_zero(critical_dir1)) {
        return result;
    }

    VectorS2 critical_point1 = critical_dir1.normalized();
    VectorS2 critical_point2 = critical_dir2.normalized();

    if (arc.contains(critical_point1)) {
        result = ThetaInterval::hull(result, ThetaInterval(theta(critical_point1), 0));
    }

    if (arc.contains(critical_point2)) {
        result = ThetaInterval::hull(result, ThetaInterval(theta(critical_point2), 0));
    }

    return result;
}

inline Interval BoundingBoxCalculator::z_interval_from_endpoints(const Arc& arc) {
    return Interval::hull(arc.source().z(), arc.target().z());
}

inline double BoundingBoxCalculator::theta(const VectorS2& point) {
    double t = std::atan2(point.y(), point.x());
    return t < 0.0 ? t + 2.0 * M_PI : t;
}

inline bool BoundingBoxCalculator::is_zero(const VectorS2& v) {
    return v.squaredNorm() < GEOMETRIC_EPSILON;
}

inline bool BoundingBoxCalculator::is_zero(double v) {
    return std::abs(v) < GEOMETRIC_EPSILON;
}

inline VectorS2 BoundingBoxCalculator::project_to_xy_plane(const VectorS2& v) {
    return VectorS2(v.x(), v.y(), 0.0);
}

inline VectorS2 BoundingBoxCalculator::perpendicular_in_xy_plane(const VectorS2& v) {
    return VectorS2(v.y(), -v.x(), 0.0);
}

inline bool BoundingBoxCalculator::parallel(const VectorS2& v1, const VectorS2& v2) {
    double v1_squared_length = v1.squaredNorm();
    double v2_squared_length = v2.squaredNorm();

    if (v1_squared_length < GEOMETRIC_EPSILON || v2_squared_length < GEOMETRIC_EPSILON) {
        return true;
    }

    VectorS2 cross = v1.cross(v2);
    return cross.squaredNorm() / (v1_squared_length * v2_squared_length) < GEOMETRIC_EPSILON;
}

} // namespace geometry_art::geometry::spherical::polygon

namespace geometry_art::geometry::spherical {
using PolygonBoundingBoxCalculator = polygon::BoundingBoxCalculator;
}

namespace geometry_art {
using PolygonBoundingBoxCalculator = geometry::spherical::polygon::BoundingBoxCalculator;
}

#endif //GEOMETRY_ART_GEOMETRY_SPHERICAL_POLYGON_BOUNDING_BOX_CALCULATOR_HPP_
