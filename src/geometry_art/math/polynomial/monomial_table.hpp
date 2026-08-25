#ifndef GEOMETRY_ART_MATH_POLYNOMIAL_MONOMIAL_TABLE_HPP_
#define GEOMETRY_ART_MATH_POLYNOMIAL_MONOMIAL_TABLE_HPP_

#include "multi_index.hpp"
#include <CGAL/assertions.h>
#include <cstddef>
#include <vector>

namespace geometry_art::math::polynomial {

class MonomialTable {
 public:
    explicit MonomialTable(int max_degree);

    [[nodiscard]] int max_degree() const { return _max_degree; }
    [[nodiscard]] size_t size() const { return _values.size(); }
    [[nodiscard]] const std::vector<double>& values() const { return _values; }

    [[nodiscard]] double at(const MultiIndex& index) const;
    [[nodiscard]] double at(int x, int y, int z) const { return at(MultiIndex{x, y, z}); }
    void set(const MultiIndex& index, double value);
    void add(const MultiIndex& index, double value);

    [[nodiscard]] bool operator==(const MonomialTable& other) const = default;

 private:
    int _max_degree;
    std::vector<double> _values;
};

inline MonomialTable::MonomialTable(int max_degree) :
    _max_degree(max_degree),
    _values(MultiIndex::count_up_to(max_degree), 0.0) {
    CGAL_precondition(max_degree >= 0);
}

inline double MonomialTable::at(const MultiIndex& index) const {
    CGAL_precondition(index.degree() <= _max_degree);
    return _values[MultiIndex::ordinal(index)];
}

inline void MonomialTable::set(const MultiIndex& index, double value) {
    CGAL_precondition(index.degree() <= _max_degree);
    _values[MultiIndex::ordinal(index)] = value;
}

inline void MonomialTable::add(const MultiIndex& index, double value) {
    CGAL_precondition(index.degree() <= _max_degree);
    _values[MultiIndex::ordinal(index)] += value;
}

} // namespace geometry_art::math::polynomial

#endif //GEOMETRY_ART_MATH_POLYNOMIAL_MONOMIAL_TABLE_HPP_
