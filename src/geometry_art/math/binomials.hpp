#ifndef GEOMETRY_ART_MATH_BINOMIALS_HPP_
#define GEOMETRY_ART_MATH_BINOMIALS_HPP_

#include "polynomial/maximum_degree.hpp"
#include <CGAL/assertions.h>
#include <vector>

namespace geometry_art::math {

// Pascal's triangle indexed [n][k], covering every degree the caller asked
// for. Expanding a boundary parametrisation into monomials needs these for
// every degree at once, so they are read from a table rather than evaluated
// per term; a row depends on nothing but its own index, so one table built
// to the maximum degree serves every caller.
[[nodiscard]] inline const std::vector<std::vector<double>>& binomials([[maybe_unused]] int max_degree) {
    CGAL_precondition(max_degree >= 0 && max_degree <= polynomial::MAXIMUM_DEGREE);

    static const std::vector<std::vector<double>> table = [] {
        int bound = polynomial::MAXIMUM_DEGREE;
        std::vector<std::vector<double>> result(bound + 1, std::vector<double>(bound + 1, 0.0));

        for (int n = 0; n <= bound; ++n) {
            result[n][0] = 1.0;

            for (int k = 1; k <= n; ++k) {
                result[n][k] = result[n - 1][k - 1] + (k <= n - 1 ? result[n - 1][k] : 0.0);
            }
        }

        return result;
    }();

    return table;
}

} // namespace geometry_art::math

#endif //GEOMETRY_ART_MATH_BINOMIALS_HPP_
