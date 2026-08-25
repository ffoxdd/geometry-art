#ifndef GEOMETRY_ART_MATH_POWERS_HPP_
#define GEOMETRY_ART_MATH_POWERS_HPP_

#include <vector>

namespace geometry_art::math {

// Successive powers of a base, indexed by exponent. Expanding a
// parametrisation into monomials needs every exponent up to a degree, and
// each follows from the last by one multiplication, so a general power is
// never called.
[[nodiscard]] inline std::vector<double> powers(double base, int max_degree) {
    std::vector<double> result(max_degree + 1, 1.0);

    for (int degree = 1; degree <= max_degree; ++degree) {
        result[degree] = result[degree - 1] * base;
    }

    return result;
}

} // namespace geometry_art::math

#endif //GEOMETRY_ART_MATH_POWERS_HPP_
