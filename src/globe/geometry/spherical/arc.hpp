#ifndef GLOBEART_SRC_GLOBE_GEOMETRY_SPHERICAL_ARC_HPP_
#define GLOBEART_SRC_GLOBE_GEOMETRY_SPHERICAL_ARC_HPP_

#include "../../math/binomials.hpp"
#include "helpers.hpp"
#include "../../types.hpp"
#include "../../math/polynomial/moments.hpp"
#include "../../math/polynomial/multi_index.hpp"
#include <cstdint>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include <optional>
#include <vector>

namespace globe::geometry::spherical {

using globe::VectorS2;
using globe::GEOMETRIC_EPSILON;
using globe::math::polynomial::Moments;
using globe::math::polynomial::MultiIndex;

class Arc {
 public:
    Arc(const VectorS2& source, const VectorS2& target, const VectorS2& normal);
    Arc(const VectorS2& source, const VectorS2& target);

    [[nodiscard]] const VectorS2& source() const { return _source; }
    [[nodiscard]] const VectorS2& target() const { return _target; }
    [[nodiscard]] const VectorS2& normal() const { return _normal; }
    [[nodiscard]] VectorS2 tangent_at_source() const;

    [[nodiscard]] double length() const;
    [[nodiscard]] VectorS2 interpolate(double t) const;
    [[nodiscard]] Arc subarc(const VectorS2& point) const;
    [[nodiscard]] bool contains(const VectorS2& point) const;
    [[nodiscard]] std::optional<Arc> clipped_by(const VectorS2& half_space_normal) const;
    [[nodiscard]] VectorS2 crossing_with(const VectorS2& great_circle_normal) const;

    [[nodiscard]] Moments moments(int max_degree) const;
    [[nodiscard]] VectorS2 first_moment() const;
    [[nodiscard]] Eigen::Matrix3d second_moment() const;

 private:
    VectorS2 _source;
    VectorS2 _target;
    VectorS2 _normal;

    [[nodiscard]] double angle_to(const VectorS2& point) const;

    [[nodiscard]] static VectorS2 normal_for(const VectorS2& source, const VectorS2& target);
    [[nodiscard]] static VectorS2 perpendicular_to(const VectorS2& u);
    [[nodiscard]] static std::vector<std::vector<double>> trigonometric_integrals(double theta, int max_degree);
    [[nodiscard]] static std::vector<double> powers(double base, int max_degree);
};

inline Arc::Arc(
    const VectorS2& source,
    const VectorS2& target,
    const VectorS2& normal
) : _source(source), _target(target), _normal(normal) {
}

inline Arc::Arc(const VectorS2& source, const VectorS2& target) :
    _source(source), _target(target), _normal(normal_for(source, target)) {
}

inline VectorS2 Arc::tangent_at_source() const {
    return _normal.cross(_source);
}

inline double Arc::length() const {
    return angle_to(_target);
}

inline VectorS2 Arc::interpolate(double t) const {
    double angle = t * length();
    return (std::cos(angle) * _source + std::sin(angle) * tangent_at_source()).normalized();
}

inline Arc Arc::subarc(const VectorS2& point) const {
    return Arc(_source, point, _normal);
}

inline bool Arc::contains(const VectorS2& point) const {
    if (std::abs(_normal.dot(point)) > GEOMETRIC_EPSILON) {
        return false;
    }

    double angle = angle_to(point);
    return angle <= length() + GEOMETRIC_EPSILON || angle >= TWO_PI - GEOMETRIC_EPSILON;
}

inline std::optional<Arc> Arc::clipped_by(const VectorS2& half_space_normal) const {
    bool source_inside = half_space_normal.dot(_source) >= -GEOMETRIC_EPSILON;
    bool target_inside = half_space_normal.dot(_target) >= -GEOMETRIC_EPSILON;

    if (source_inside && target_inside) {
        return *this;
    }

    if (!source_inside && !target_inside) {
        return std::nullopt;
    }

    VectorS2 crossing = crossing_with(half_space_normal);
    return source_inside ? Arc(_source, crossing, _normal) : Arc(crossing, _target, _normal);
}

inline VectorS2 Arc::crossing_with(const VectorS2& great_circle_normal) const {
    VectorS2 candidate = _normal.cross(great_circle_normal).normalized();
    double theta = length();

    auto excess = [&](const VectorS2& point) {
        double angle = angle_to(point);
        return std::min(std::max(0.0, angle - theta), TWO_PI - angle);
    };

    return excess(candidate) <= excess(-candidate) ? candidate : VectorS2(-candidate);
}

inline double Arc::angle_to(const VectorS2& point) const {
    double angle = std::atan2(tangent_at_source().dot(point), _source.dot(point));
    return angle < 0.0 ? angle + TWO_PI : angle;
}

inline Moments Arc::moments(int max_degree) const {
    VectorS2 u = _source;
    VectorS2 n = tangent_at_source();
    auto integrals = trigonometric_integrals(length(), max_degree);
    auto binomial = globe::math::binomials(max_degree);
    std::vector<std::vector<double>> u_powers{powers(u.x(), max_degree), powers(u.y(), max_degree), powers(u.z(), max_degree)};
    std::vector<std::vector<double>> n_powers{powers(n.x(), max_degree), powers(n.y(), max_degree), powers(n.z(), max_degree)};

    Moments result(max_degree);

    for (const MultiIndex& alpha : MultiIndex::all_up_to(max_degree)) {
        double total = 0.0;

        for (int beta_x = 0; beta_x <= alpha.x; ++beta_x) {
            for (int beta_y = 0; beta_y <= alpha.y; ++beta_y) {
                for (int beta_z = 0; beta_z <= alpha.z; ++beta_z) {
                    int cosine_power = beta_x + beta_y + beta_z;
                    int sine_power = alpha.degree() - cosine_power;

                    double weight =
                        binomial[alpha.x][beta_x] * u_powers[0][beta_x] * n_powers[0][alpha.x - beta_x] *
                        binomial[alpha.y][beta_y] * u_powers[1][beta_y] * n_powers[1][alpha.y - beta_y] *
                        binomial[alpha.z][beta_z] * u_powers[2][beta_z] * n_powers[2][alpha.z - beta_z];

                    total += weight * integrals[cosine_power][sine_power];
                }
            }
        }

        result.set(alpha, total);
    }

    return result;
}

inline VectorS2 Arc::first_moment() const {
    Moments moments = this->moments(1);
    return VectorS2(moments.at(1, 0, 0), moments.at(0, 1, 0), moments.at(0, 0, 1));
}

inline Eigen::Matrix3d Arc::second_moment() const {
    Moments moments = this->moments(2);
    Eigen::Matrix3d result;

    for (int row = 0; row < 3; ++row) {
        for (int column = 0; column < 3; ++column) {
            result(row, column) = moments.at(MultiIndex::unit(row).raised(column));
        }
    }

    return result;
}

inline VectorS2 Arc::normal_for(const VectorS2& source, const VectorS2& target) {
    VectorS2 cross = source.cross(target);
    double norm = cross.norm();

    if (norm < GEOMETRIC_EPSILON) {
        return perpendicular_to(source);
    }

    return cross / norm;
}

inline VectorS2 Arc::perpendicular_to(const VectorS2& u) {
    VectorS2 candidate = (std::abs(u.z()) < 0.9) ? VectorS2(0, 0, 1) : VectorS2(1, 0, 0);
    return (candidate - u.dot(candidate) * u).normalized();
}

inline std::vector<std::vector<double>> Arc::trigonometric_integrals(double theta, int max_degree) {
    double sine = std::sin(theta);
    double cosine = std::cos(theta);
    auto sine_powers = powers(sine, max_degree + 1);
    auto cosine_powers = powers(cosine, max_degree + 1);

    std::vector<std::vector<double>> integrals(
        max_degree + 1,
        std::vector<double>(max_degree + 1, 0.0)
    );

    integrals[0][0] = theta;

    for (int total = 1; total <= max_degree; ++total) {
        for (int p = 0; p <= total; ++p) {
            int q = total - p;

            if (p == 1 && q == 0) {
                integrals[p][q] = sine;
            } else if (p == 0 && q == 1) {
                integrals[p][q] = 1.0 - cosine;
            } else if (p == 1 && q == 1) {
                integrals[p][q] = sine * sine / 2.0;
            } else if (p >= 2) {
                integrals[p][q] =
                    (cosine_powers[p - 1] * sine_powers[q + 1] + (p - 1) * integrals[p - 2][q]) / total;
            } else {
                integrals[p][q] =
                    (-cosine_powers[p + 1] * sine_powers[q - 1] + (q - 1) * integrals[p][q - 2]) / total;
            }
        }
    }

    return integrals;
}

inline std::vector<double> Arc::powers(double base, int max_degree) {
    std::vector<double> result(max_degree + 1, 1.0);

    for (int degree = 1; degree <= max_degree; ++degree) {
        result[degree] = result[degree - 1] * base;
    }

    return result;
}


} // namespace globe::geometry::spherical

namespace globe {
using Arc = geometry::spherical::Arc;
}

#endif //GLOBEART_SRC_GLOBE_GEOMETRY_SPHERICAL_ARC_HPP_
