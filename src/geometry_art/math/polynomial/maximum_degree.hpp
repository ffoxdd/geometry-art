#ifndef GEOMETRY_ART_MATH_POLYNOMIAL_MAXIMUM_DEGREE_HPP_
#define GEOMETRY_ART_MATH_POLYNOMIAL_MAXIMUM_DEGREE_HPP_

namespace geometry_art::math::polynomial {

// The highest degree the moment machinery supports. Tables that depend on
// nothing but the degree are built once up to this bound and shared, so the
// bound is what makes that sharing finite.
constexpr int MAXIMUM_DEGREE = 32;

} // namespace geometry_art::math::polynomial

#endif //GEOMETRY_ART_MATH_POLYNOMIAL_MAXIMUM_DEGREE_HPP_
