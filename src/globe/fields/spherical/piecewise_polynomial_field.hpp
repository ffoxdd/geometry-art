#ifndef GLOBEART_SRC_GLOBE_FIELDS_SPHERICAL_PIECEWISE_POLYNOMIAL_FIELD_HPP_
#define GLOBEART_SRC_GLOBE_FIELDS_SPHERICAL_PIECEWISE_POLYNOMIAL_FIELD_HPP_

#include "field.hpp"
#include "polynomial_field.hpp"
#include "region_integrals.hpp"
#include "../scalar/field.hpp"
#include "../../types.hpp"
#include "../../cgal/types.hpp"
#include "../../geometry/spherical/arc.hpp"
#include "../../geometry/spherical/helpers.hpp"
#include "../../geometry/spherical/indexed_kd_tree.hpp"
#include "../../geometry/spherical/polygon/polygon.hpp"
#include "../../geometry/spherical/triangle_mesh.hpp"
#include "../../math/polynomial/moments.hpp"
#include "../../math/polynomial/multi_index.hpp"
#include "../../math/polynomial/polynomial.hpp"
#include <CGAL/assertions.h>
#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

namespace globe::fields::spherical {

using geometry::spherical::IndexedFuzzySphere;
using geometry::spherical::IndexedKDTree;
using geometry::spherical::IndexedPointMap;
using geometry::spherical::IndexedSearchTraits;
using geometry::spherical::TriangleMesh;
using geometry::spherical::distance;
using globe::math::polynomial::Moments;
using globe::math::polynomial::MultiIndex;
using globe::math::polynomial::Polynomial;

class PiecewisePolynomialField {
 public:
    PiecewisePolynomialField(TriangleMesh mesh, std::vector<Polynomial> triangle_polynomials);

    template<scalar::Field ScalarFieldType>
    [[nodiscard]] static PiecewisePolynomialField sample(TriangleMesh mesh, int degree, ScalarFieldType& scalar_field);

    [[nodiscard]] const TriangleMesh& mesh() const { return _mesh; }
    [[nodiscard]] int degree() const { return _degree; }

    [[nodiscard]] double value(const VectorS2& point) const;
    [[nodiscard]] RegionIntegrals integrals(const Polygon& polygon) const;
    [[nodiscard]] RegionIntegrals integrals(const Arc& arc) const;
    [[nodiscard]] RegionIntegrals integrals(const Polygon& polygon, const std::vector<Moments>& arc_moments) const;
    [[nodiscard]] RegionIntegrals integrals(const Arc& arc, const Moments& arc_moments) const;
    [[nodiscard]] double total_mass() const { return _total_mass; }

    [[nodiscard]] static std::vector<VectorS2> lagrange_nodes(const std::array<VectorS2, 3>& corners, int degree);
    [[nodiscard]] static Polynomial interpolate(const std::array<VectorS2, 3>& corners, int degree, const std::vector<double>& node_values);

 private:
    struct Triangle {
        std::array<VectorS2, 3> inward_normals;
        PolynomialField field;
        VectorS2 center;
        double angular_radius;
    };

    struct SpatialIndex {
        std::vector<cgal::Point3> centers;
        IndexedPointMap<std::vector<cgal::Point3>> point_map;
        IndexedSearchTraits traits;
        IndexedKDTree tree;

        explicit SpatialIndex(std::vector<cgal::Point3> triangle_centers);
    };

    TriangleMesh _mesh;
    int _degree;
    std::vector<Triangle> _triangles;
    std::shared_ptr<SpatialIndex> _index;
    double _max_angular_radius;
    double _total_mass;

    [[nodiscard]] std::vector<size_t> candidates(const VectorS2& center, double angular_radius) const;
    [[nodiscard]] static std::array<VectorS2, 3> corners_of(const TriangleMesh& mesh, const std::array<size_t, 3>& indices);
    [[nodiscard]] static Triangle build_triangle(const std::array<VectorS2, 3>& corners, Polynomial polynomial);
    [[nodiscard]] static Polygon triangle_polygon(const std::array<VectorS2, 3>& corners);
    [[nodiscard]] static double monomial(const VectorS2& point, const MultiIndex& index);
};

inline PiecewisePolynomialField::SpatialIndex::SpatialIndex(std::vector<cgal::Point3> triangle_centers) :
    centers(std::move(triangle_centers)),
    point_map(centers),
    traits(point_map),
    tree(
        boost::counting_iterator<size_t>(0),
        boost::counting_iterator<size_t>(centers.size()),
        IndexedKDTree::Splitter(),
        traits
    ) {
    tree.build();
}

inline PiecewisePolynomialField::PiecewisePolynomialField(TriangleMesh mesh, std::vector<Polynomial> triangle_polynomials) :
    _mesh(std::move(mesh)),
    _degree(triangle_polynomials.empty() ? 0 : triangle_polynomials.front().max_degree()),
    _max_angular_radius(0.0),
    _total_mass(0.0) {
    CGAL_precondition(triangle_polynomials.size() == _mesh.triangles.size());

    std::vector<cgal::Point3> centers;
    _triangles.reserve(_mesh.triangles.size());
    centers.reserve(_mesh.triangles.size());

    for (size_t i = 0; i < _mesh.triangles.size(); ++i) {
        std::array<VectorS2, 3> corners = corners_of(_mesh, _mesh.triangles[i]);
        _triangles.push_back(build_triangle(corners, std::move(triangle_polynomials[i])));
        centers.push_back(cgal::to_point(_triangles.back().center));
        _max_angular_radius = std::max(_max_angular_radius, _triangles.back().angular_radius);
        _total_mass += _triangles.back().field.integrals(triangle_polygon(corners)).mass;
    }

    _index = std::make_shared<SpatialIndex>(std::move(centers));
}

template<scalar::Field ScalarFieldType>
PiecewisePolynomialField PiecewisePolynomialField::sample(TriangleMesh mesh, int degree, ScalarFieldType& scalar_field) {
    std::vector<Polynomial> polynomials;
    polynomials.reserve(mesh.triangles.size());

    for (const auto& indices : mesh.triangles) {
        std::array<VectorS2, 3> corners = corners_of(mesh, indices);
        std::vector<double> values;

        for (const VectorS2& node : lagrange_nodes(corners, degree)) {
            values.push_back(scalar_field.value(node));
        }

        polynomials.push_back(interpolate(corners, degree, values));
    }

    return PiecewisePolynomialField(std::move(mesh), std::move(polynomials));
}

inline std::vector<VectorS2> PiecewisePolynomialField::lagrange_nodes(const std::array<VectorS2, 3>& corners, int degree) {
    std::vector<VectorS2> nodes;

    for (int i = degree; i >= 0; --i) {
        for (int j = degree - i; j >= 0; --j) {
            int k = degree - i - j;
            nodes.push_back((i * corners[0] + j * corners[1] + k * corners[2]).normalized());
        }
    }

    return nodes;
}

inline Polynomial PiecewisePolynomialField::interpolate(
    const std::array<VectorS2, 3>& corners,
    int degree,
    const std::vector<double>& node_values
) {
    std::vector<VectorS2> nodes = lagrange_nodes(corners, degree);
    std::vector<MultiIndex> basis = MultiIndex::all_of_degree(degree);
    CGAL_precondition(node_values.size() == nodes.size());
    CGAL_precondition(basis.size() == nodes.size());

    Eigen::MatrixXd design(nodes.size(), basis.size());
    Eigen::VectorXd values(nodes.size());

    for (size_t row = 0; row < nodes.size(); ++row) {
        values[static_cast<Eigen::Index>(row)] = node_values[row];
        for (size_t column = 0; column < basis.size(); ++column) {
            design(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column)) = monomial(nodes[row], basis[column]);
        }
    }

    Eigen::VectorXd coefficients = design.fullPivLu().solve(values);
    Polynomial polynomial(degree);

    for (size_t i = 0; i < basis.size(); ++i) {
        polynomial.set_coefficient(basis[i], coefficients[static_cast<Eigen::Index>(i)]);
    }

    return polynomial;
}

inline double PiecewisePolynomialField::value(const VectorS2& point) const {
    for (size_t index : candidates(point, 0.0)) {
        const Triangle& triangle = _triangles[index];
        bool inside = true;

        for (const VectorS2& normal : triangle.inward_normals) {
            inside = inside && normal.dot(point) >= -GEOMETRIC_EPSILON;
        }

        if (inside) {
            return triangle.field.value(point);
        }
    }

    return 0.0;
}

inline RegionIntegrals PiecewisePolynomialField::integrals(const Polygon& polygon) const {
    VectorS2 center = polygon.centroid();
    double angular_radius = 2.0 * std::asin(std::min(1.0, polygon.bounding_sphere_radius() / 2.0));
    RegionIntegrals total{0.0, Vector3::Zero()};

    for (size_t index : candidates(center, angular_radius)) {
        const Triangle& triangle = _triangles[index];
        std::optional<Polygon> piece = polygon;

        for (const VectorS2& normal : triangle.inward_normals) {
            if (!piece) {
                break;
            }
            piece = piece->clipped_by(normal);
        }

        if (piece) {
            RegionIntegrals contribution = triangle.field.integrals(*piece);
            total.mass += contribution.mass;
            total.first_moment += contribution.first_moment;
        }
    }

    return total;
}

inline RegionIntegrals PiecewisePolynomialField::integrals(const Arc& arc) const {
    VectorS2 center = arc.interpolate(0.5);
    double angular_radius = arc.length() / 2.0;
    RegionIntegrals total{0.0, Vector3::Zero()};

    for (size_t index : candidates(center, angular_radius)) {
        const Triangle& triangle = _triangles[index];
        std::optional<Arc> piece = arc;

        for (const VectorS2& normal : triangle.inward_normals) {
            if (!piece) {
                break;
            }
            piece = piece->clipped_by(normal);
        }

        if (piece) {
            RegionIntegrals contribution = triangle.field.integrals(*piece);
            total.mass += contribution.mass;
            total.first_moment += contribution.first_moment;
        }
    }

    return total;
}

inline RegionIntegrals PiecewisePolynomialField::integrals(
    const Polygon& polygon,
    [[maybe_unused]] const std::vector<Moments>& arc_moments
) const {
    return integrals(polygon);
}

inline RegionIntegrals PiecewisePolynomialField::integrals(
    const Arc& arc,
    [[maybe_unused]] const Moments& arc_moments
) const {
    return integrals(arc);
}

inline std::vector<size_t> PiecewisePolynomialField::candidates(const VectorS2& center, double angular_radius) const {
    double search_angle = std::min(M_PI, angular_radius + _max_angular_radius);
    double chord_radius = 2.0 * std::sin(search_angle / 2.0);
    std::vector<size_t> result;

    IndexedFuzzySphere query(cgal::to_point(center), chord_radius, GEOMETRIC_EPSILON, _index->traits);
    _index->tree.search(std::back_inserter(result), query);

    return result;
}

inline std::array<VectorS2, 3> PiecewisePolynomialField::corners_of(const TriangleMesh& mesh, const std::array<size_t, 3>& indices) {
    return {mesh.vertices[indices[0]], mesh.vertices[indices[1]], mesh.vertices[indices[2]]};
}

inline PiecewisePolynomialField::Triangle PiecewisePolynomialField::build_triangle(
    const std::array<VectorS2, 3>& corners,
    Polynomial polynomial
) {
    VectorS2 center = (corners[0] + corners[1] + corners[2]).normalized();
    double angular_radius = 0.0;
    for (const VectorS2& corner : corners) {
        angular_radius = std::max(angular_radius, distance(center, corner));
    }

    return Triangle{
        {
            corners[0].cross(corners[1]).normalized(),
            corners[1].cross(corners[2]).normalized(),
            corners[2].cross(corners[0]).normalized()
        },
        PolynomialField(std::move(polynomial)),
        center,
        angular_radius
    };
}

inline Polygon PiecewisePolynomialField::triangle_polygon(const std::array<VectorS2, 3>& corners) {
    return Polygon(std::vector<Arc>{
        Arc(corners[0], corners[1]),
        Arc(corners[1], corners[2]),
        Arc(corners[2], corners[0])
    });
}

inline double PiecewisePolynomialField::monomial(const VectorS2& point, const MultiIndex& index) {
    return std::pow(point.x(), index.x) * std::pow(point.y(), index.y) * std::pow(point.z(), index.z);
}

static_assert(Field<PiecewisePolynomialField>);

} // namespace globe::fields::spherical

#endif //GLOBEART_SRC_GLOBE_FIELDS_SPHERICAL_PIECEWISE_POLYNOMIAL_FIELD_HPP_
