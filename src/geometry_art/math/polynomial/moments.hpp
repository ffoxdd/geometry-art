#ifndef GEOMETRY_ART_MATH_POLYNOMIAL_MOMENTS_HPP_
#define GEOMETRY_ART_MATH_POLYNOMIAL_MOMENTS_HPP_

#include "monomial_table.hpp"
#include "multi_index.hpp"
#include <cmath>

namespace geometry_art::math::polynomial {

class Moments {
 public:
    explicit Moments(int max_degree);

    [[nodiscard]] int max_degree() const { return _table.max_degree(); }
    [[nodiscard]] const MonomialTable& table() const { return _table; }

    [[nodiscard]] double at(const MultiIndex& index) const { return _table.at(index); }
    [[nodiscard]] double at(int x, int y, int z) const { return _table.at(x, y, z); }
    void set(const MultiIndex& index, double value) { _table.set(index, value); }
    void add(const MultiIndex& index, double value) { _table.add(index, value); }

    [[nodiscard]] static Moments unit_sphere(int max_degree);

 private:
    MonomialTable _table;

    [[nodiscard]] static double unit_sphere_monomial_integral(const MultiIndex& index);
    [[nodiscard]] static double double_factorial(int n);
};

inline Moments::Moments(int max_degree) :
    _table(max_degree) {
}

inline Moments Moments::unit_sphere(int max_degree) {
    Moments moments(max_degree);

    for (const MultiIndex& index : MultiIndex::all_up_to(max_degree)) {
        moments.set(index, unit_sphere_monomial_integral(index));
    }

    return moments;
}

inline double Moments::unit_sphere_monomial_integral(const MultiIndex& index) {
    if (index.x % 2 != 0 || index.y % 2 != 0 || index.z % 2 != 0) {
        return 0.0;
    }

    double numerator =
        double_factorial(index.x - 1) *
        double_factorial(index.y - 1) *
        double_factorial(index.z - 1);

    return 4.0 * M_PI * numerator / double_factorial(index.degree() + 1);
}

inline double Moments::double_factorial(int n) {
    double result = 1.0;

    for (int k = n; k > 1; k -= 2) {
        result *= k;
    }

    return result;
}

} // namespace geometry_art::math::polynomial

#endif //GEOMETRY_ART_MATH_POLYNOMIAL_MOMENTS_HPP_
