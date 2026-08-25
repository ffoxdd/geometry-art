#ifndef GEOMETRY_ART_CGAL_TYPES_HPP_
#define GEOMETRY_ART_CGAL_TYPES_HPP_

#include "../types.hpp"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>

namespace geometry_art::cgal {

template<typename T>
concept VectorXYZ = geometry_art::VectorXYZ<T> && requires(const T& v) {
    ::CGAL::to_double(v.x());
    ::CGAL::to_double(v.y());
    ::CGAL::to_double(v.z());
};

using Kernel = ::CGAL::Exact_predicates_inexact_constructions_kernel;
using Vector3 = Kernel::Vector_3;
using Point3 = Kernel::Point_3;
using SearchTraits = ::CGAL::Search_traits_3<Kernel>;
using KDTree = ::CGAL::Kd_tree<SearchTraits>;
using FuzzySphere = ::CGAL::Fuzzy_sphere<SearchTraits>;

template<geometry_art::VectorXYZ VectorType>
inline Point3 to_point(const VectorType& v) {
    return Point3(v.x(), v.y(), v.z());
}

template<geometry_art::VectorXYZ VectorType>
inline Vector3 to_vector(const VectorType& v) {
    return Vector3(v.x(), v.y(), v.z());
}

}

namespace geometry_art {

    template<cgal::VectorXYZ CGALVectorType>
    inline Vector3 to_vector3(const CGALVectorType& v) {
        return Vector3(
            ::CGAL::to_double(v.x()),
            ::CGAL::to_double(v.y()),
            ::CGAL::to_double(v.z())
        );
    }

    template<cgal::VectorXYZ CGALVectorType>
    inline VectorS2 to_vector_s2(const CGALVectorType& v) {
        return VectorS2(
            ::CGAL::to_double(v.x()),
            ::CGAL::to_double(v.y()),
            ::CGAL::to_double(v.z())
        );
    }

    template<cgal::VectorXYZ CGALVectorType>
    inline Eigen::Vector3d to_eigen(const CGALVectorType& v) {
        return Eigen::Vector3d(
            ::CGAL::to_double(v.x()),
            ::CGAL::to_double(v.y()),
            ::CGAL::to_double(v.z())
        );
    }

}

#endif //GEOMETRY_ART_CGAL_TYPES_HPP_
