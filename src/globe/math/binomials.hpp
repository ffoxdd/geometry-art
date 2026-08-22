#ifndef GLOBEART_SRC_GLOBE_MATH_BINOMIALS_HPP_
#define GLOBEART_SRC_GLOBE_MATH_BINOMIALS_HPP_

#include <vector>

namespace globe::math {

// Pascal's triangle up to max_degree, indexed [n][k]. Expanding a boundary
// parametrisation into monomials needs these for every degree at once, so
// they are built as a table rather than evaluated per term.
[[nodiscard]] inline std::vector<std::vector<double>> binomials(int max_degree) {
    std::vector<std::vector<double>> result(max_degree + 1, std::vector<double>(max_degree + 1, 0.0));

    for (int n = 0; n <= max_degree; ++n) {
        result[n][0] = 1.0;

        for (int k = 1; k <= n; ++k) {
            result[n][k] = result[n - 1][k - 1] + (k <= n - 1 ? result[n - 1][k] : 0.0);
        }
    }

    return result;
}

} // namespace globe::math

#endif //GLOBEART_SRC_GLOBE_MATH_BINOMIALS_HPP_
