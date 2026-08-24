#ifndef GLOBEART_SRC_GLOBE_FIELDS_FLAT_PIECEWISE_POLYNOMIAL_FIELD_HPP_
#define GLOBEART_SRC_GLOBE_FIELDS_FLAT_PIECEWISE_POLYNOMIAL_FIELD_HPP_

#include "scalar_field.hpp"
#include "../region_integrals.hpp"
#include "../../types.hpp"
#include "../../geometry/planar/polygon.hpp"
#include "../../geometry/planar/segment.hpp"
#include "../../math/polynomial/moments.hpp"
#include "../../math/polynomial/polynomial.hpp"
#include <CGAL/assertions.h>
#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

namespace globe::fields::flat {

using geometry::planar::Polygon;
using geometry::planar::Segment;
using globe::math::polynomial::Moments;
using globe::math::polynomial::Polynomial;

// A density held as one quadratic per triangle of a regular periodic grid
// over the rectangle of periods: each grid square splits into two
// triangles, and each triangle carries the quadratic through its corner and
// edge-midpoint samples. Shared nodes are sampled at one canonical
// position, so the field is continuous, and continuous across the seams.
//
// A query region lives in some cell's chart and may protrude past the
// rectangle, so integration walks the period translates its bounding box
// touches, clips in mesh coordinates, and translates the resulting moments
// back into the chart.
class PiecewisePolynomialField {
 public:
    template<ScalarField ScalarFieldType>
    [[nodiscard]] static PiecewisePolynomialField sample(
        double width,
        double height,
        int columns,
        int rows,
        ScalarFieldType& scalar_field
    );

    [[nodiscard]] double value(const Vector2& point) const;
    [[nodiscard]] int degree() const { return DEGREE; }
    [[nodiscard]] double total_mass() const { return _total_mass; }
    [[nodiscard]] double width() const { return _width; }
    [[nodiscard]] double height() const { return _height; }
    [[nodiscard]] size_t piece_count() const { return _pieces.size(); }
    [[nodiscard]] double lowest_sampled_value() const { return _lowest_sampled_value; }

    [[nodiscard]] RegionIntegrals integrals(const Polygon& polygon) const;
    [[nodiscard]] RegionIntegrals integrals(const Segment& segment) const;
    [[nodiscard]] double squared_norm_moment(const Polygon& polygon) const;
    [[nodiscard]] Matrix3 second_moment(const Segment& segment) const;
    [[nodiscard]] Vector3 gradient_masses(const Segment& segment) const;
    [[nodiscard]] Matrix3 gradient_first_moments(const Segment& segment) const;
    [[nodiscard]] std::array<Matrix3, 3> gradient_second_moments(const Segment& segment) const;

 private:
    static constexpr int DEGREE = 2;
    static constexpr int NODE_COUNT = 6;

    double _width;
    double _height;
    int _columns;
    int _rows;
    double _total_mass;
    double _lowest_sampled_value;
    std::vector<Polynomial> _pieces;
    std::vector<std::array<Polynomial, 2>> _pieces_times_coordinate;
    std::vector<Polynomial> _pieces_times_squared_norm;
    std::vector<std::array<Polynomial, 3>> _pieces_times_coordinate_pair;
    std::vector<std::array<Polynomial, 2>> _piece_derivatives;
    std::vector<std::array<std::array<Polynomial, 2>, 2>> _piece_derivatives_times_coordinate;
    std::vector<std::array<std::array<Polynomial, 3>, 2>> _piece_derivatives_times_coordinate_pair;

    PiecewisePolynomialField(double width, double height, int columns, int rows);

    [[nodiscard]] std::array<Vector2, 3> triangle_corners(size_t piece) const;
    [[nodiscard]] std::optional<Polygon> clipped_to_triangle(const Polygon& polygon, size_t piece) const;
    [[nodiscard]] std::optional<Segment> clipped_to_triangle(const Segment& segment, size_t piece) const;

    template<typename RegionType, typename Accumulate>
    void for_each_overlap(const RegionType& region, const Accumulate& accumulate) const;

    [[nodiscard]] static Polynomial fit_piece(
        const std::array<Vector2, 3>& corners,
        const std::vector<double>& node_values
    );
    [[nodiscard]] static std::array<double, 2> bounding_interval(
        const std::vector<Vector2>& points,
        int axis
    );
};

template<ScalarField ScalarFieldType>
PiecewisePolynomialField PiecewisePolynomialField::sample(
    double width,
    double height,
    int columns,
    int rows,
    ScalarFieldType& scalar_field
) {
    CGAL_precondition(width > 0.0 && height > 0.0 && columns > 0 && rows > 0);

    PiecewisePolynomialField field(width, height, columns, rows);
    field._pieces.reserve(static_cast<size_t>(columns) * rows * 2);

    auto canonical_value = [&](const Vector2& point) {
        double x = std::fmod(point.x(), width);
        double y = std::fmod(point.y(), height);

        if (x < 0.0) { x += width; }
        if (y < 0.0) { y += height; }

        return scalar_field.value(Vector2(x, y));
    };

    for (size_t piece = 0; piece < static_cast<size_t>(columns) * rows * 2; ++piece) {
        std::array<Vector2, 3> corners = field.triangle_corners(piece);
        std::vector<double> node_values;
        node_values.reserve(NODE_COUNT);

        for (int corner = 0; corner < 3; ++corner) {
            node_values.push_back(canonical_value(corners[corner]));
        }

        for (int corner = 0; corner < 3; ++corner) {
            node_values.push_back(canonical_value(0.5 * (corners[corner] + corners[(corner + 1) % 3])));
        }

        field._lowest_sampled_value = std::min(
            field._lowest_sampled_value,
            *std::min_element(node_values.begin(), node_values.end())
        );

        field._pieces.push_back(fit_piece(corners, node_values));
    }

    field._pieces_times_coordinate.reserve(field._pieces.size());
    field._pieces_times_squared_norm.reserve(field._pieces.size());
    field._pieces_times_coordinate_pair.reserve(field._pieces.size());

    for (size_t piece = 0; piece < field._pieces.size(); ++piece) {
        const Polynomial& density = field._pieces[piece];
        Polynomial times_x = density.times_coordinate(0);
        Polynomial times_y = density.times_coordinate(1);

        field._pieces_times_coordinate.push_back({times_x, times_y});
        field._pieces_times_squared_norm.push_back(
            times_x.times_coordinate(0).plus(times_y.times_coordinate(1))
        );
        field._pieces_times_coordinate_pair.push_back({
            times_x.times_coordinate(0),
            times_x.times_coordinate(1),
            times_y.times_coordinate(1)
        });

        auto coordinate_products = [](const Polynomial& derivative) {
            return std::array<Polynomial, 2>{derivative.times_coordinate(0), derivative.times_coordinate(1)};
        };
        auto pair_products = [](const std::array<Polynomial, 2>& products) {
            return std::array<Polynomial, 3>{
                products[0].times_coordinate(0),
                products[0].times_coordinate(1),
                products[1].times_coordinate(1)
            };
        };

        std::array<Polynomial, 2> derivatives{
            density.partial_derivative(0),
            density.partial_derivative(1)
        };
        std::array<std::array<Polynomial, 2>, 2> derivative_coordinates{
            coordinate_products(derivatives[0]),
            coordinate_products(derivatives[1])
        };

        field._piece_derivatives_times_coordinate_pair.push_back({
            pair_products(derivative_coordinates[0]),
            pair_products(derivative_coordinates[1])
        });
        field._piece_derivatives.push_back(std::move(derivatives));
        field._piece_derivatives_times_coordinate.push_back(std::move(derivative_coordinates));

        std::array<Vector2, 3> corners = field.triangle_corners(piece);
        Moments moments = Polygon(std::vector<Vector2>{corners[0], corners[1], corners[2]}).moments(DEGREE);
        field._total_mass += density.integrate(moments);
    }

    return field;
}

inline PiecewisePolynomialField::PiecewisePolynomialField(double width, double height, int columns, int rows) :
    _width(width),
    _height(height),
    _columns(columns),
    _rows(rows),
    _total_mass(0.0),
    _lowest_sampled_value(std::numeric_limits<double>::infinity()) {
}

// Piece order: two per grid square, the lower-left triangle first, columns
// fastest. The lower-left triangle holds the corner of its square; the
// upper-right holds the opposite one.
inline std::array<Vector2, 3> PiecewisePolynomialField::triangle_corners(size_t piece) const {
    size_t square = piece / 2;
    int column = static_cast<int>(square % static_cast<size_t>(_columns));
    int row = static_cast<int>(square / static_cast<size_t>(_columns));
    double cell_width = _width / _columns;
    double cell_height = _height / _rows;

    Vector2 low(column * cell_width, row * cell_height);
    Vector2 right = low + Vector2(cell_width, 0.0);
    Vector2 up = low + Vector2(0.0, cell_height);
    Vector2 high = low + Vector2(cell_width, cell_height);

    if (piece % 2 == 0) {
        return {low, right, up};
    }

    return {right, high, up};
}

inline double PiecewisePolynomialField::value(const Vector2& point) const {
    double x = std::fmod(point.x(), _width);
    double y = std::fmod(point.y(), _height);

    if (x < 0.0) { x += _width; }
    if (y < 0.0) { y += _height; }

    double cell_width = _width / _columns;
    double cell_height = _height / _rows;
    int column = std::min(static_cast<int>(x / cell_width), _columns - 1);
    int row = std::min(static_cast<int>(y / cell_height), _rows - 1);
    double fraction_x = x / cell_width - column;
    double fraction_y = y / cell_height - row;

    size_t piece = 2 * (static_cast<size_t>(row) * _columns + column) + (fraction_x + fraction_y <= 1.0 ? 0 : 1);

    return _pieces[piece].value(Vector3(x, y, 0.0));
}

// Walks every (piece, period translate) pair whose triangle can meet the
// region, clips in mesh coordinates, and hands the clipped region, the
// piece and the translate to the accumulator.
template<typename RegionType, typename Accumulate>
void PiecewisePolynomialField::for_each_overlap(const RegionType& region, const Accumulate& accumulate) const {
    std::vector<Vector2> points;

    if constexpr (std::is_same_v<RegionType, Polygon>) {
        points = region.vertices();
    } else {
        points = {region.source(), region.target()};
    }

    std::array<double, 2> x_range = bounding_interval(points, 0);
    std::array<double, 2> y_range = bounding_interval(points, 1);

    double cell_width = _width / _columns;
    double cell_height = _height / _rows;

    for (int tile_x = static_cast<int>(std::floor(x_range[0] / _width));
         tile_x <= static_cast<int>(std::floor(x_range[1] / _width)); ++tile_x) {
        for (int tile_y = static_cast<int>(std::floor(y_range[0] / _height));
             tile_y <= static_cast<int>(std::floor(y_range[1] / _height)); ++tile_y) {
            Vector2 offset(tile_x * _width, tile_y * _height);

            int column_low = std::max(0, static_cast<int>(std::floor((x_range[0] - offset.x()) / cell_width)));
            int column_high = std::min(_columns - 1, static_cast<int>(std::floor((x_range[1] - offset.x()) / cell_width)));
            int row_low = std::max(0, static_cast<int>(std::floor((y_range[0] - offset.y()) / cell_height)));
            int row_high = std::min(_rows - 1, static_cast<int>(std::floor((y_range[1] - offset.y()) / cell_height)));

            for (int row = row_low; row <= row_high; ++row) {
                for (int column = column_low; column <= column_high; ++column) {
                    for (int half = 0; half < 2; ++half) {
                        size_t piece = 2 * (static_cast<size_t>(row) * _columns + column) + half;
                        accumulate(piece, offset);
                    }
                }
            }
        }
    }
}

inline RegionIntegrals PiecewisePolynomialField::integrals(const Polygon& polygon) const {
    RegionIntegrals result{0.0, Vector3::Zero()};

    for_each_overlap(polygon, [&](size_t piece, const Vector2& offset) {
        std::vector<Vector2> shifted;
        shifted.reserve(polygon.size());

        for (const Vector2& vertex : polygon.vertices()) {
            shifted.push_back(vertex - offset);
        }

        std::optional<Polygon> clipped = clipped_to_triangle(Polygon(std::move(shifted)), piece);

        if (!clipped.has_value()) {
            return;
        }

        Moments moments = clipped->moments(DEGREE + 1);
        double mass = _pieces[piece].integrate(moments);
        Vector3 first(
            _pieces_times_coordinate[piece][0].integrate(moments),
            _pieces_times_coordinate[piece][1].integrate(moments),
            0.0
        );

        result.mass += mass;
        result.first_moment += first + Vector3(offset.x(), offset.y(), 0.0) * mass;
    });

    return result;
}

inline RegionIntegrals PiecewisePolynomialField::integrals(const Segment& segment) const {
    RegionIntegrals result{0.0, Vector3::Zero()};

    for_each_overlap(segment, [&](size_t piece, const Vector2& offset) {
        std::optional<Segment> clipped = clipped_to_triangle(
            Segment(segment.source() - offset, segment.target() - offset),
            piece
        );

        if (!clipped.has_value()) {
            return;
        }

        Moments moments = clipped->moments(DEGREE + 1);
        double mass = _pieces[piece].integrate(moments);
        Vector3 first(
            _pieces_times_coordinate[piece][0].integrate(moments),
            _pieces_times_coordinate[piece][1].integrate(moments),
            0.0
        );

        result.mass += mass;
        result.first_moment += first + Vector3(offset.x(), offset.y(), 0.0) * mass;
    });

    return result;
}

inline double PiecewisePolynomialField::squared_norm_moment(const Polygon& polygon) const {
    double result = 0.0;

    for_each_overlap(polygon, [&](size_t piece, const Vector2& offset) {
        std::vector<Vector2> shifted;
        shifted.reserve(polygon.size());

        for (const Vector2& vertex : polygon.vertices()) {
            shifted.push_back(vertex - offset);
        }

        std::optional<Polygon> clipped = clipped_to_triangle(Polygon(std::move(shifted)), piece);

        if (!clipped.has_value()) {
            return;
        }

        Moments moments = clipped->moments(DEGREE + 2);
        double mass = _pieces[piece].integrate(moments);
        double squared = _pieces_times_squared_norm[piece].integrate(moments);
        Vector2 first(
            _pieces_times_coordinate[piece][0].integrate(moments),
            _pieces_times_coordinate[piece][1].integrate(moments)
        );

        result += squared + 2.0 * offset.dot(first) + offset.squaredNorm() * mass;
    });

    return result;
}

inline Matrix3 PiecewisePolynomialField::second_moment(const Segment& segment) const {
    Matrix3 result = Matrix3::Zero();

    for_each_overlap(segment, [&](size_t piece, const Vector2& offset) {
        std::optional<Segment> clipped = clipped_to_triangle(
            Segment(segment.source() - offset, segment.target() - offset),
            piece
        );

        if (!clipped.has_value()) {
            return;
        }

        Moments moments = clipped->moments(DEGREE + 2);
        double mass = _pieces[piece].integrate(moments);
        Vector3 first(
            _pieces_times_coordinate[piece][0].integrate(moments),
            _pieces_times_coordinate[piece][1].integrate(moments),
            0.0
        );

        Matrix3 mesh = Matrix3::Zero();
        mesh(0, 0) = _pieces_times_coordinate_pair[piece][0].integrate(moments);
        mesh(0, 1) = _pieces_times_coordinate_pair[piece][1].integrate(moments);
        mesh(1, 0) = mesh(0, 1);
        mesh(1, 1) = _pieces_times_coordinate_pair[piece][2].integrate(moments);

        Vector3 shift(offset.x(), offset.y(), 0.0);
        result += mesh + shift * first.transpose() + first * shift.transpose() +
            mass * shift * shift.transpose();
    });

    return result;
}

inline Vector3 PiecewisePolynomialField::gradient_masses(const Segment& segment) const {
    Vector3 result = Vector3::Zero();

    for_each_overlap(segment, [&](size_t piece, const Vector2& offset) {
        std::optional<Segment> clipped = clipped_to_triangle(
            Segment(segment.source() - offset, segment.target() - offset),
            piece
        );

        if (!clipped.has_value()) {
            return;
        }

        Moments moments = clipped->moments(DEGREE + 1);

        for (int axis = 0; axis < 2; ++axis) {
            result[axis] += _piece_derivatives[piece][axis].integrate(moments);
        }
    });

    return result;
}

inline Matrix3 PiecewisePolynomialField::gradient_first_moments(const Segment& segment) const {
    Matrix3 result = Matrix3::Zero();

    for_each_overlap(segment, [&](size_t piece, const Vector2& offset) {
        std::optional<Segment> clipped = clipped_to_triangle(
            Segment(segment.source() - offset, segment.target() - offset),
            piece
        );

        if (!clipped.has_value()) {
            return;
        }

        Moments moments = clipped->moments(DEGREE + 1);
        Vector3 shift(offset.x(), offset.y(), 0.0);

        for (int axis = 0; axis < 2; ++axis) {
            double mass = _piece_derivatives[piece][axis].integrate(moments);
            Vector3 first(
                _piece_derivatives_times_coordinate[piece][axis][0].integrate(moments),
                _piece_derivatives_times_coordinate[piece][axis][1].integrate(moments),
                0.0
            );

            result.col(axis) += first + shift * mass;
        }
    });

    return result;
}

inline std::array<Matrix3, 3> PiecewisePolynomialField::gradient_second_moments(const Segment& segment) const {
    std::array<Matrix3, 3> result{Matrix3::Zero(), Matrix3::Zero(), Matrix3::Zero()};

    for_each_overlap(segment, [&](size_t piece, const Vector2& offset) {
        std::optional<Segment> clipped = clipped_to_triangle(
            Segment(segment.source() - offset, segment.target() - offset),
            piece
        );

        if (!clipped.has_value()) {
            return;
        }

        Moments moments = clipped->moments(DEGREE + 2);
        Vector3 shift(offset.x(), offset.y(), 0.0);

        for (int axis = 0; axis < 2; ++axis) {
            double mass = _piece_derivatives[piece][axis].integrate(moments);
            Vector3 first(
                _piece_derivatives_times_coordinate[piece][axis][0].integrate(moments),
                _piece_derivatives_times_coordinate[piece][axis][1].integrate(moments),
                0.0
            );

            Matrix3 mesh = Matrix3::Zero();
            mesh(0, 0) = _piece_derivatives_times_coordinate_pair[piece][axis][0].integrate(moments);
            mesh(0, 1) = _piece_derivatives_times_coordinate_pair[piece][axis][1].integrate(moments);
            mesh(1, 0) = mesh(0, 1);
            mesh(1, 1) = _piece_derivatives_times_coordinate_pair[piece][axis][2].integrate(moments);

            result[axis] += mesh + shift * first.transpose() + first * shift.transpose() +
                mass * shift * shift.transpose();
        }
    });

    return result;
}

inline std::optional<Polygon> PiecewisePolynomialField::clipped_to_triangle(
    const Polygon& polygon,
    size_t piece
) const {
    std::array<Vector2, 3> corners = triangle_corners(piece);
    std::optional<Polygon> clipped = polygon;

    for (int side = 0; side < 3 && clipped.has_value(); ++side) {
        Vector2 edge = corners[(side + 1) % 3] - corners[side];
        clipped = clipped->clipped_by(Vector2(-edge.y(), edge.x()), corners[side]);
    }

    return clipped;
}

inline std::optional<Segment> PiecewisePolynomialField::clipped_to_triangle(
    const Segment& segment,
    size_t piece
) const {
    std::array<Vector2, 3> corners = triangle_corners(piece);
    double low = 0.0;
    double high = 1.0;

    for (int side = 0; side < 3; ++side) {
        Vector2 edge = corners[(side + 1) % 3] - corners[side];
        Vector2 inward(-edge.y(), edge.x());
        double at_source = inward.dot(segment.source() - corners[side]);
        double at_target = inward.dot(segment.target() - corners[side]);
        double span = at_target - at_source;

        if (std::abs(span) < GEOMETRIC_EPSILON * std::max(1.0, std::abs(at_source))) {
            if (at_source < 0.0) {
                return std::nullopt;
            }

            continue;
        }

        double crossing = -at_source / span;

        if (span > 0.0) {
            low = std::max(low, crossing);
        } else {
            high = std::min(high, crossing);
        }
    }

    if (low >= high) {
        return std::nullopt;
    }

    return Segment(segment.interpolate(low), segment.interpolate(high));
}

inline Polynomial PiecewisePolynomialField::fit_piece(
    const std::array<Vector2, 3>& corners,
    const std::vector<double>& node_values
) {
    std::vector<Vector2> nodes;
    nodes.reserve(NODE_COUNT);

    for (int corner = 0; corner < 3; ++corner) {
        nodes.push_back(corners[corner]);
    }

    for (int corner = 0; corner < 3; ++corner) {
        nodes.push_back(0.5 * (corners[corner] + corners[(corner + 1) % 3]));
    }

    Eigen::Matrix<double, NODE_COUNT, NODE_COUNT> design;
    Eigen::Matrix<double, NODE_COUNT, 1> values;

    for (int node = 0; node < NODE_COUNT; ++node) {
        double x = nodes[node].x();
        double y = nodes[node].y();
        design.row(node) << 1.0, x, y, x * x, x * y, y * y;
        values[node] = node_values[node];
    }

    Eigen::Matrix<double, NODE_COUNT, 1> coefficients = design.colPivHouseholderQr().solve(values);

    Polynomial piece(DEGREE);
    piece.add_coefficient({0, 0, 0}, coefficients[0]);
    piece.add_coefficient({1, 0, 0}, coefficients[1]);
    piece.add_coefficient({0, 1, 0}, coefficients[2]);
    piece.add_coefficient({2, 0, 0}, coefficients[3]);
    piece.add_coefficient({1, 1, 0}, coefficients[4]);
    piece.add_coefficient({0, 2, 0}, coefficients[5]);

    return piece;
}

inline std::array<double, 2> PiecewisePolynomialField::bounding_interval(
    const std::vector<Vector2>& points,
    int axis
) {
    double low = std::numeric_limits<double>::infinity();
    double high = -std::numeric_limits<double>::infinity();

    for (const Vector2& point : points) {
        low = std::min(low, point[axis]);
        high = std::max(high, point[axis]);
    }

    return {low, high};
}

} // namespace globe::fields::flat

#endif //GLOBEART_SRC_GLOBE_FIELDS_FLAT_PIECEWISE_POLYNOMIAL_FIELD_HPP_
