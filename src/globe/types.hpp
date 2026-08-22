#ifndef GLOBEART_SRC_GLOBE_TYPES_H_
#define GLOBEART_SRC_GLOBE_TYPES_H_

#include <cstdint>
#include <Eigen/Core>

namespace globe {

template<typename T>
concept VectorXYZ = requires(const T &value) { value.x(); value.y(); value.z(); };

// Canonical type for 3D Cartesian points.
using Vector3 = Eigen::Vector3d;
using Matrix3 = Eigen::Matrix3d;

// Canonical type for points on the unit sphere.
// Represents a normalized 3D vector (direction from origin to point on S2).
using VectorS2 = Eigen::Vector3d;

constexpr double TWO_PI = 2.0 * M_PI;
constexpr double GEOMETRIC_EPSILON = 1e-10;

template<VectorXYZ VectorType>
inline Vector3 to_vector3(const VectorType& p) {
    return Vector3(p.x(), p.y(), p.z());
}

template<VectorXYZ VectorType>
inline VectorS2 to_vector_s2(const VectorType& p) {
    return VectorS2(p.x(), p.y(), p.z());
}

template<VectorXYZ VectorType>
inline Eigen::Vector3d to_eigen(const VectorType& p) {
    return Eigen::Vector3d(p.x(), p.y(), p.z());
}

}

#endif //GLOBEART_SRC_GLOBE_TYPES_H_
