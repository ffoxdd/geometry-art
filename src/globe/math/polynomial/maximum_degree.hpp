#ifndef GLOBEART_SRC_GLOBE_MATH_POLYNOMIAL_MAXIMUM_DEGREE_HPP_
#define GLOBEART_SRC_GLOBE_MATH_POLYNOMIAL_MAXIMUM_DEGREE_HPP_

namespace globe::math::polynomial {

// The highest degree the moment machinery supports. Tables that depend on
// nothing but the degree are built once up to this bound and shared, so the
// bound is what makes that sharing finite.
constexpr int MAXIMUM_DEGREE = 32;

} // namespace globe::math::polynomial

#endif //GLOBEART_SRC_GLOBE_MATH_POLYNOMIAL_MAXIMUM_DEGREE_HPP_
