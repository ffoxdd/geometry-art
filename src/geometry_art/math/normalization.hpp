#ifndef GEOMETRY_ART_MATH_NORMALIZATION_HPP_
#define GEOMETRY_ART_MATH_NORMALIZATION_HPP_

#include "../types.hpp"

namespace geometry_art::math {

// Derivatives of f(y) = h(y / ||y||) expressed in the derivatives of h at
// the normalised point. Sites are optimised as ambient vectors and used as
// directions, so every gradient and Hessian the optimiser sees passes
// through this map.
class Normalization {
 public:
    explicit Normalization(const Vector3& point);

    [[nodiscard]] const Vector3& site() const { return _site; }
    [[nodiscard]] double norm() const { return _norm; }

    [[nodiscard]] Vector3 gradient(const Vector3& site_gradient) const;
    [[nodiscard]] Matrix3 hessian(const Vector3& site_gradient, const Matrix3& site_hessian) const;
    [[nodiscard]] Matrix3 tangential_hessian(const Vector3& site_gradient, const Matrix3& site_hessian) const;
    [[nodiscard]] Matrix3 mixed_hessian(const Matrix3& site_hessian, const Normalization& column) const;

 private:
    Vector3 _site;
    double _norm;
    Matrix3 _projection;
};

inline Normalization::Normalization(const Vector3& point) :
    _site(point.normalized()),
    _norm(point.norm()),
    _projection(Matrix3::Identity() - _site * _site.transpose()) {
}

inline Vector3 Normalization::gradient(const Vector3& site_gradient) const {
    return _projection * site_gradient / _norm;
}

inline Matrix3 Normalization::hessian(const Vector3& site_gradient, const Matrix3& site_hessian) const {
    Vector3 projected = _projection * site_gradient;
    double radial = _site.dot(site_gradient);

    return (
        _projection * site_hessian * _projection -
        radial * _projection -
        projected * _site.transpose() -
        _site * projected.transpose()
    ) / (_norm * _norm);
}

// Scaling a site changes nothing, so the radial direction carries no
// curvature and would read to a solver as a direction worth exploring.
// Restricting to the tangent plane removes it.
inline Matrix3 Normalization::tangential_hessian(
    const Vector3& site_gradient,
    const Matrix3& site_hessian
) const {
    return _projection * hessian(site_gradient, site_hessian) * _projection;
}

inline Matrix3 Normalization::mixed_hessian(const Matrix3& site_hessian, const Normalization& column) const {
    return _projection * site_hessian * column._projection / (_norm * column._norm);
}

} // namespace geometry_art::math

namespace geometry_art {
using Normalization = math::Normalization;
}

#endif //GEOMETRY_ART_MATH_NORMALIZATION_HPP_
