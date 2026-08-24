#ifndef GLOBEART_SRC_GLOBE_FIELDS_SPHERICAL_POWELL_SABIN_PROJECTION_HPP_
#define GLOBEART_SRC_GLOBE_FIELDS_SPHERICAL_POWELL_SABIN_PROJECTION_HPP_

#include "piecewise_polynomial_field.hpp"
#include "../scalar/field.hpp"
#include "../../types.hpp"
#include "../../geometry/spherical/barycentric.hpp"
#include "../../geometry/spherical/bezier_triangle.hpp"
#include "../../geometry/spherical/powell_sabin_refinement.hpp"
#include "../../geometry/spherical/triangle_mesh.hpp"
#include "../../math/polynomial/multi_index.hpp"
#include "../../math/polynomial/polynomial.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

namespace globe::fields::spherical {

using geometry::spherical::Barycentric;
using geometry::spherical::BezierTriangle;
using geometry::spherical::PowellSabinRefinement;
using geometry::spherical::TriangleMesh;
using globe::math::polynomial::MultiIndex;
using globe::math::polynomial::Polynomial;

// Builds a C1 field from any callable as the least-squares projection of the
// callable onto the spline space, fitted over dense samples across the whole
// sphere in one sparse symmetric solve.
//
// The pieces are quadratics on the Powell-Sabin refinement, which is the
// coarsest split on which a C1 quadratic spline exists. Three consequences
// matter here:
//
//   - the density's gradient is continuous, so the constraint curvature the
//     optimizer needs is continuous too, and the second-order method has the
//     smoothness its convergence rests on;
//   - the pieces are bounded by their Bezier coefficients, so a lower bound
//     on the density is read off the coefficients rather than sampled for,
//     and a fit that dips below zero is impossible to miss;
//   - the field is the best the space can do in the mean-square sense, so
//     structure finer than the mesh lands in the residual instead of being
//     aliased into the result, and the local averages the tessellation reads
//     are as faithful as the mesh allows.
//
// The space is parametrized by one vector per mesh vertex. A piece is a
// quadratic homogeneous in space, so by Euler's relation its value and its
// tangential gradient at a vertex are together one gradient vector
// `h = 2 f V + g`, and the spline is linear in these: three unknowns per
// vertex, and a projection that is one sparse positive-definite system.
//
// Every piece is still a polynomial on a great-circle triangle, so the
// integration path is the one already in use.
class PowellSabinProjection {
 public:
    struct Result {
        PiecewisePolynomialField field;

        // The least Bezier coefficient over every piece. The field is at
        // least this everywhere, so a positive value certifies a positive
        // density without sampling.
        double lowest_coefficient;

        // The smallest gradient damping any vertex needed. One means the
        // projection was left as solved; less than one means it dipped below
        // zero somewhere, and some vertex gradients were reduced to keep the
        // density positive.
        double least_damping;

        // The root mean square of what the space could not represent,
        // against the sphere's area: the measured representation error, and
        // the number a requested accuracy would be checked against.
        double root_mean_square_residual;
    };

    // Each sub-triangle of the split is sampled on a grid this many rows
    // deep, one sample per grid triangle, so the fit sees the field between
    // the vertices and not only at them.
    explicit PowellSabinProjection(int sample_side = DEFAULT_SAMPLE_SIDE);

    template<scalar::Field ScalarFieldType>
    [[nodiscard]] Result project(const TriangleMesh& mesh, ScalarFieldType& scalar_field) const;

 private:
    // The domain points of the six quadratics, deduplicated: the seven
    // vertices of the split, then the midpoint of each of its twelve edges.
    static constexpr int DOMAIN_COUNT = 19;
    static constexpr int CORNER_DATA_COUNT = 9;
    static constexpr int VERTEX_CONDITION_COUNT = 12;
    static constexpr int SMOOTHNESS_CONDITION_COUNT = 12;
    static constexpr int DEGREE = 2;
    static constexpr int DEFAULT_SAMPLE_SIDE = 4;
    static constexpr int MAXIMUM_DAMPING_ROUNDS = 60;
    static constexpr double DAMPING_FACTOR = 0.5;

    using CoefficientMap = Eigen::Matrix<double, DOMAIN_COUNT, CORNER_DATA_COUNT>;

    static constexpr std::array<std::array<int, 3>, 12> EDGE_MIDPOINTS = {{
        {PowellSabinRefinement::CORNER_0, PowellSabinRefinement::SPLIT_01, 7},
        {PowellSabinRefinement::SPLIT_01, PowellSabinRefinement::CORNER_1, 8},
        {PowellSabinRefinement::CORNER_1, PowellSabinRefinement::SPLIT_12, 9},
        {PowellSabinRefinement::SPLIT_12, PowellSabinRefinement::CORNER_2, 10},
        {PowellSabinRefinement::CORNER_2, PowellSabinRefinement::SPLIT_20, 11},
        {PowellSabinRefinement::SPLIT_20, PowellSabinRefinement::CORNER_0, 12},
        {PowellSabinRefinement::CORNER_0, PowellSabinRefinement::INTERIOR, 13},
        {PowellSabinRefinement::CORNER_1, PowellSabinRefinement::INTERIOR, 14},
        {PowellSabinRefinement::CORNER_2, PowellSabinRefinement::INTERIOR, 15},
        {PowellSabinRefinement::SPLIT_01, PowellSabinRefinement::INTERIOR, 16},
        {PowellSabinRefinement::SPLIT_12, PowellSabinRefinement::INTERIOR, 17},
        {PowellSabinRefinement::SPLIT_20, PowellSabinRefinement::INTERIOR, 18}
    }};

    int _sample_side;

    struct Projection {
        Eigen::VectorXd solution;
        double root_mean_square_residual;
    };

    template<scalar::Field ScalarFieldType>
    [[nodiscard]] Projection solve_projection(
        const PowellSabinRefinement& refinement,
        const std::vector<CoefficientMap>& maps,
        ScalarFieldType& scalar_field
    ) const;

    [[nodiscard]] static CoefficientMap coefficients_from_corner_data(const PowellSabinRefinement::Cell& cell);

    [[nodiscard]] static Eigen::Matrix<double, Eigen::Dynamic, DOMAIN_COUNT> bernstein_rows(
        const Barycentric& barycentric,
        const std::array<int, 3>& sub_triangle,
        const std::vector<VectorS2>& points
    );

    [[nodiscard]] static std::vector<Vector3> barycentric_samples(int side);
    [[nodiscard]] static double spherical_area(const Vector3& a, const Vector3& b, const Vector3& c);

    [[nodiscard]] static Eigen::Matrix<double, CORNER_DATA_COUNT, 1> corner_data(
        const PowellSabinRefinement::Cell& cell,
        const std::vector<double>& values,
        const std::vector<Vector3>& gradients,
        const std::vector<double>& damping
    );

    [[nodiscard]] static std::array<Vector3, DOMAIN_COUNT> domain_points(const PowellSabinRefinement::Cell& cell);
    [[nodiscard]] static int domain_index(int first, int second);
    [[nodiscard]] static int domain_index_for(const std::array<int, 3>& corner_order, const MultiIndex& index);
    [[nodiscard]] static bool contains(const std::array<int, 3>& values, int value);
};

inline PowellSabinProjection::PowellSabinProjection(int sample_side) :
    _sample_side(sample_side) {
    CGAL_precondition(sample_side >= 1);
}

template<scalar::Field ScalarFieldType>
PowellSabinProjection::Result PowellSabinProjection::project(
    const TriangleMesh& mesh,
    ScalarFieldType& scalar_field
) const {
    PowellSabinRefinement refinement(mesh);
    const std::vector<PowellSabinRefinement::Cell>& cells = refinement.cells();

    std::vector<CoefficientMap> maps;
    maps.reserve(cells.size());

    for (const PowellSabinRefinement::Cell& cell : cells) {
        maps.push_back(coefficients_from_corner_data(cell));
    }

    Projection projection = solve_projection(refinement, maps, scalar_field);
    const Eigen::VectorXd& solution = projection.solution;

    // By Euler's relation the solved vector at a vertex carries the value in
    // its radial part and the tangential gradient in the rest.
    std::vector<double> values(mesh.vertices.size());
    std::vector<Vector3> gradients(mesh.vertices.size());

    for (size_t index = 0; index < mesh.vertices.size(); ++index) {
        Vector3 vertex_data = solution.segment<3>(3 * static_cast<Eigen::Index>(index));
        double radial = vertex_data.dot(mesh.vertices[index]);
        values[index] = 0.5 * radial;
        gradients[index] = vertex_data - radial * mesh.vertices[index];
    }

    std::vector<double> damping(mesh.vertices.size(), 1.0);
    std::vector<Eigen::VectorXd> coefficients(cells.size());
    std::vector<bool> stale(cells.size(), true);

    // A projection is optimal in the mean square, not bounded by its target,
    // so it can dip below zero where the mesh is far too coarse for the
    // field. Damping is per vertex rather than per piece, so every piece
    // meeting at a vertex still sees the same data and the field stays C1
    // while it is brought back into range.
    for (int round = 0; round <= MAXIMUM_DAMPING_ROUNDS; ++round) {
        for (size_t index = 0; index < cells.size(); ++index) {
            if (stale[index]) {
                coefficients[index] = maps[index] * corner_data(cells[index], values, gradients, damping);
                stale[index] = false;
            }
        }

        std::vector<bool> overshoots(mesh.vertices.size(), false);
        bool any = false;

        for (size_t index = 0; index < cells.size(); ++index) {
            if (coefficients[index].minCoeff() < 0.0) {
                any = true;

                for (size_t corner : cells[index].corner_vertices) {
                    overshoots[corner] = true;
                }
            }
        }

        if (!any || round == MAXIMUM_DAMPING_ROUNDS) {
            break;
        }

        for (size_t index = 0; index < mesh.vertices.size(); ++index) {
            if (overshoots[index]) {
                damping[index] *= DAMPING_FACTOR;
            }
        }

        for (size_t index = 0; index < cells.size(); ++index) {
            for (size_t corner : cells[index].corner_vertices) {
                stale[index] = stale[index] || overshoots[corner];
            }
        }
    }

    std::vector<Polynomial> polynomials;
    polynomials.reserve(refinement.mesh().triangles.size());
    double lowest = std::numeric_limits<double>::infinity();

    for (size_t index = 0; index < cells.size(); ++index) {
        for (const std::array<int, 3>& sub_triangle : PowellSabinRefinement::SUB_TRIANGLES) {
            BezierTriangle piece(
                {
                    cells[index].points[sub_triangle[0]],
                    cells[index].points[sub_triangle[1]],
                    cells[index].points[sub_triangle[2]]
                },
                DEGREE
            );

            for (const MultiIndex& multi_index : BezierTriangle::indices(DEGREE)) {
                piece.set_coefficient(multi_index, coefficients[index][domain_index_for(sub_triangle, multi_index)]);
            }

            lowest = std::min(lowest, piece.lowest_coefficient());
            polynomials.push_back(piece.to_polynomial());
        }
    }

    return Result{
        PiecewisePolynomialField(refinement.mesh(), std::move(polynomials)),
        lowest,
        *std::min_element(damping.begin(), damping.end()),
        projection.root_mean_square_residual
    };
}

// The normal equations of a least-squares fit over samples spread across
// every sub-triangle, weighted by area so the fit approximates the L2
// projection. Each sample touches only the three vertices of its cell, so
// the system is sparse, symmetric and positive definite.
template<scalar::Field ScalarFieldType>
typename PowellSabinProjection::Projection PowellSabinProjection::solve_projection(
    const PowellSabinRefinement& refinement,
    const std::vector<CoefficientMap>& maps,
    ScalarFieldType& scalar_field
) const {
    const std::vector<PowellSabinRefinement::Cell>& cells = refinement.cells();
    size_t vertex_count = 0;

    for (const PowellSabinRefinement::Cell& cell : cells) {
        for (size_t corner : cell.corner_vertices) {
            vertex_count = std::max(vertex_count, corner + 1);
        }
    }

    Eigen::Index unknown_count = 3 * static_cast<Eigen::Index>(vertex_count);
    std::vector<Eigen::Triplet<double>> triplets;
    Eigen::VectorXd right_hand_side = Eigen::VectorXd::Zero(unknown_count);
    std::vector<Vector3> samples = barycentric_samples(_sample_side);
    double squared_target = 0.0;
    double total_weight = 0.0;

    for (size_t index = 0; index < cells.size(); ++index) {
        const PowellSabinRefinement::Cell& cell = cells[index];
        Eigen::Matrix<double, CORNER_DATA_COUNT, CORNER_DATA_COUNT> normal_block;
        Eigen::Matrix<double, CORNER_DATA_COUNT, 1> data_block;
        normal_block.setZero();
        data_block.setZero();

        for (const std::array<int, 3>& sub_triangle : PowellSabinRefinement::SUB_TRIANGLES) {
            Barycentric barycentric({
                cell.points[sub_triangle[0]],
                cell.points[sub_triangle[1]],
                cell.points[sub_triangle[2]]
            });

            std::vector<VectorS2> points;
            points.reserve(samples.size());

            for (const Vector3& sample : samples) {
                points.push_back(VectorS2(barycentric.combination(sample)).normalized());
            }

            Eigen::Matrix<double, Eigen::Dynamic, DOMAIN_COUNT> bernstein =
                bernstein_rows(barycentric, sub_triangle, points);

            Eigen::Matrix<double, Eigen::Dynamic, CORNER_DATA_COUNT> design = bernstein * maps[index];
            Eigen::VectorXd sampled(static_cast<Eigen::Index>(points.size()));

            for (size_t sample = 0; sample < points.size(); ++sample) {
                sampled[static_cast<Eigen::Index>(sample)] = scalar_field.value(points[sample]);
            }

            double weight = spherical_area(
                cell.points[sub_triangle[0]],
                cell.points[sub_triangle[1]],
                cell.points[sub_triangle[2]]
            ) / static_cast<double>(samples.size());

            normal_block += weight * design.transpose() * design;
            data_block += weight * design.transpose() * sampled;
            squared_target += weight * sampled.squaredNorm();
            total_weight += weight * static_cast<double>(points.size());
        }

        for (int row = 0; row < CORNER_DATA_COUNT; ++row) {
            Eigen::Index global_row = 3 * static_cast<Eigen::Index>(cell.corner_vertices[row / 3]) + row % 3;
            right_hand_side[global_row] += data_block[row];

            for (int column = 0; column < CORNER_DATA_COUNT; ++column) {
                Eigen::Index global_column =
                    3 * static_cast<Eigen::Index>(cell.corner_vertices[column / 3]) + column % 3;
                triplets.emplace_back(global_row, global_column, normal_block(row, column));
            }
        }
    }

    Eigen::SparseMatrix<double> normal(unknown_count, unknown_count);
    normal.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(normal);
    CGAL_postcondition(solver.info() == Eigen::Success);
    Eigen::VectorXd solution = solver.solve(right_hand_side);

    // At the optimum the squared residual telescopes to what the samples
    // carry minus what the solution reproduces, so the representation error
    // costs nothing extra to know.
    double squared_residual = std::max(0.0, squared_target - solution.dot(right_hand_side));

    return Projection{std::move(solution), std::sqrt(squared_residual / total_weight)};
}

// One small dense solve per cell, mapping the corner vectors to the
// coefficients at the nineteen domain points. The rows say what the corners
// pin down and that the pieces meet smoothly inside the cell; smoothness
// across a cell's boundary comes from the refinement's geometry, not from
// these rows. Since every condition is linear in the corner vectors, the
// whole map is one matrix, solved for all nine inputs at once.
inline PowellSabinProjection::CoefficientMap
PowellSabinProjection::coefficients_from_corner_data(const PowellSabinRefinement::Cell& cell) {
    constexpr int ROW_COUNT = VERTEX_CONDITION_COUNT + SMOOTHNESS_CONDITION_COUNT;

    // The domain points touching each corner. Setting them from the corner's
    // vector is what pins the value and the gradient there: for a quadratic
    // the coefficient beside a corner is half the vector's reading of the
    // corner it leans towards.
    constexpr std::array<std::array<int, 4>, 3> CORNER_DOMAIN_POINTS = {{
        {0, 7, 12, 13},
        {1, 8, 9, 14},
        {2, 10, 11, 15}
    }};

    std::array<Vector3, DOMAIN_COUNT> points = domain_points(cell);
    Eigen::Matrix<double, ROW_COUNT, DOMAIN_COUNT> system;
    Eigen::Matrix<double, ROW_COUNT, CORNER_DATA_COUNT> data;
    system.setZero();
    data.setZero();
    int row = 0;

    for (int corner = 0; corner < 3; ++corner) {
        for (int domain : CORNER_DOMAIN_POINTS[corner]) {
            system(row, domain) = 1.0;
            data.block<1, 3>(row, 3 * corner) = (points[domain] - 0.5 * cell.points[corner]).transpose();
            ++row;
        }
    }

    for (const std::array<int, 2>& pair : PowellSabinRefinement::ADJACENT_SUB_TRIANGLES) {
        const std::array<int, 3>& first = PowellSabinRefinement::SUB_TRIANGLES[pair[0]];
        const std::array<int, 3>& second = PowellSabinRefinement::SUB_TRIANGLES[pair[1]];

        std::array<int, 2> shared{};
        int shared_count = 0;
        int first_apex = -1;
        int second_apex = -1;

        for (int point : first) {
            if (contains(second, point)) {
                shared[shared_count++] = point;
            } else {
                first_apex = point;
            }
        }

        for (int point : second) {
            if (!contains(first, point)) {
                second_apex = point;
            }
        }

        std::array<int, 3> first_order = {first_apex, shared[0], shared[1]};
        std::array<int, 3> second_order = {second_apex, shared[0], shared[1]};

        Barycentric barycentric({cell.points[first_apex], cell.points[shared[0]], cell.points[shared[1]]});
        Vector3 apex = barycentric.coordinates(cell.points[second_apex]);

        for (const std::array<int, 2>& exponents : {std::array<int, 2>{1, 0}, std::array<int, 2>{0, 1}}) {
            int j = exponents[0];
            int k = exponents[1];

            system(row, domain_index_for(second_order, MultiIndex{1, j, k})) += 1.0;
            system(row, domain_index_for(first_order, MultiIndex{1, j, k})) -= apex[0];
            system(row, domain_index_for(first_order, MultiIndex{0, j + 1, k})) -= apex[1];
            system(row, domain_index_for(first_order, MultiIndex{0, j, k + 1})) -= apex[2];
            ++row;
        }
    }

    return system.colPivHouseholderQr().solve(data);
}

// The Bernstein basis of one sub-triangle read at each sample point,
// written against the cell's shared domain points.
inline Eigen::Matrix<double, Eigen::Dynamic, PowellSabinProjection::DOMAIN_COUNT>
PowellSabinProjection::bernstein_rows(
    const Barycentric& barycentric,
    const std::array<int, 3>& sub_triangle,
    const std::vector<VectorS2>& points
) {
    std::vector<MultiIndex> indices = BezierTriangle::indices(DEGREE);
    Eigen::Matrix<double, Eigen::Dynamic, DOMAIN_COUNT> result(
        static_cast<Eigen::Index>(points.size()),
        DOMAIN_COUNT
    );
    result.setZero();

    for (size_t sample = 0; sample < points.size(); ++sample) {
        Vector3 coordinates = barycentric.coordinates(points[sample]);

        for (const MultiIndex& index : indices) {
            double value = BezierTriangle::multinomial(index) *
                std::pow(coordinates[0], index.x) *
                std::pow(coordinates[1], index.y) *
                std::pow(coordinates[2], index.z);

            result(static_cast<Eigen::Index>(sample), domain_index_for(sub_triangle, index)) += value;
        }
    }

    return result;
}

// One sample per triangle of a regular grid over the sub-triangle, at its
// centroid, so the samples cover the piece evenly and none sits on a seam.
inline std::vector<Vector3> PowellSabinProjection::barycentric_samples(int side) {
    std::vector<Vector3> result;
    result.reserve(static_cast<size_t>(side) * static_cast<size_t>(side));
    double denominator = 3.0 * static_cast<double>(side);

    for (int i = 0; i < side; ++i) {
        for (int j = 0; i + j < side; ++j) {
            int k = side - 1 - i - j;
            result.emplace_back((3 * i + 1) / denominator, (3 * j + 1) / denominator, (3 * k + 1) / denominator);

            if (i + j < side - 1) {
                result.emplace_back((3 * i + 2) / denominator, (3 * j + 2) / denominator, (3 * k - 1) / denominator);
            }
        }
    }

    return result;
}

inline double PowellSabinProjection::spherical_area(const Vector3& a, const Vector3& b, const Vector3& c) {
    return 2.0 * std::atan2(
        std::abs(a.dot(b.cross(c))),
        1.0 + a.dot(b) + b.dot(c) + c.dot(a)
    );
}

// A piece is a quadratic that is homogeneous in space, so its radial
// derivative is twice its value: the value and the tangential gradient
// together name one vector, and every corner condition is a dot product
// against it.
inline Eigen::Matrix<double, PowellSabinProjection::CORNER_DATA_COUNT, 1> PowellSabinProjection::corner_data(
    const PowellSabinRefinement::Cell& cell,
    const std::vector<double>& values,
    const std::vector<Vector3>& gradients,
    const std::vector<double>& damping
) {
    Eigen::Matrix<double, CORNER_DATA_COUNT, 1> result;

    for (int corner = 0; corner < 3; ++corner) {
        size_t vertex = cell.corner_vertices[corner];
        result.segment<3>(3 * corner) =
            2.0 * values[vertex] * cell.points[corner] + damping[vertex] * gradients[vertex];
    }

    return result;
}

inline std::array<Vector3, PowellSabinProjection::DOMAIN_COUNT>
PowellSabinProjection::domain_points(const PowellSabinRefinement::Cell& cell) {
    std::array<Vector3, DOMAIN_COUNT> result;

    for (int point = 0; point < PowellSabinRefinement::POINT_COUNT; ++point) {
        result[point] = cell.points[point];
    }

    for (const std::array<int, 3>& midpoint : EDGE_MIDPOINTS) {
        result[midpoint[2]] = 0.5 * (cell.points[midpoint[0]] + cell.points[midpoint[1]]);
    }

    return result;
}

inline int PowellSabinProjection::domain_index(int first, int second) {
    if (first == second) {
        return first;
    }

    for (const std::array<int, 3>& midpoint : EDGE_MIDPOINTS) {
        bool forward = midpoint[0] == first && midpoint[1] == second;
        bool backward = midpoint[0] == second && midpoint[1] == first;

        if (forward || backward) {
            return midpoint[2];
        }
    }

    CGAL_error_msg("the split has no edge between these points");
    return -1;
}

inline int PowellSabinProjection::domain_index_for(const std::array<int, 3>& corner_order, const MultiIndex& index) {
    std::array<int, DEGREE> points{};
    int count = 0;

    for (int axis = 0; axis < 3; ++axis) {
        for (int repeat = 0; repeat < index[axis]; ++repeat) {
            points[count++] = corner_order[axis];
        }
    }

    return domain_index(points[0], points[1]);
}

inline bool PowellSabinProjection::contains(const std::array<int, 3>& values, int value) {
    return std::find(values.begin(), values.end(), value) != values.end();
}

} // namespace globe::fields::spherical

#endif //GLOBEART_SRC_GLOBE_FIELDS_SPHERICAL_POWELL_SABIN_PROJECTION_HPP_
