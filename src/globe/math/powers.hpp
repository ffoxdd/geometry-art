#ifndef GLOBEART_SRC_GLOBE_MATH_POWERS_HPP_
#define GLOBEART_SRC_GLOBE_MATH_POWERS_HPP_

#include <vector>

namespace globe::math {

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

} // namespace globe::math

#endif //GLOBEART_SRC_GLOBE_MATH_POWERS_HPP_
