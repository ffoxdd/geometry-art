#ifndef GLOBEART_SRC_GLOBE_FIELDS_SPHERICAL_POWELL_SABIN_INTERPOLANT_HPP_
#define GLOBEART_SRC_GLOBE_FIELDS_SPHERICAL_POWELL_SABIN_INTERPOLANT_HPP_

#include "local_quadratic_fit.hpp"
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

// Builds a C1 field from any callable by reading its value and gradient at
// the mesh vertices.
//
// The pieces are quadratics on the Powell-Sabin refinement, which is the
// coarsest split on which a C1 quadratic spline exists. Two consequences
// matter here:
//
//   - the density's gradient is continuous, so the constraint curvature the
//     optimizer needs is continuous too, and the second-order method has the
//     smoothness its convergence rests on;
//   - the pieces are bounded by their Bezier coefficients, so a lower bound
//     on the density is read off the coefficients rather than sampled for,
//     and a fit that dips below zero is impossible to miss;
//   - the vertex readings come from a fit across a neighbourhood the size
//     of the mesh rather than from samples at the vertex, so structure
//     finer than the mesh is averaged away instead of being aliased into
//     the result. That is what makes the represented field faithful in the
//     local averages the tessellation reads, at a mesh no finer than it
//     needs to be.
//
// Every piece is still a polynomial on a great-circle triangle, so the
// integration path is the one already in use.
class PowellSabinInterpolant {
 public:
    struct Result {
        PiecewisePolynomialField field;

        // The least Bezier coefficient over every piece. The field is at
        // least this everywhere, so a positive value certifies a positive
        // density without sampling.
        double lowest_coefficient;

        // The smallest gradient damping any vertex needed. One means the
        // field interpolates every sampled gradient; less than one means the
        // mesh is too coarse for how fast the density turns, and some
        // gradients were reduced to keep the density positive.
        double least_damping;
    };

    // The neighbourhood each vertex is read over, as a fraction of the mean
    // distance between neighbouring vertices. Below about a half the fit
    // stops averaging usefully; well above one it reaches past the
    // structure the mesh is able to carry.
    explicit PowellSabinInterpolant(double stencil_fraction = DEFAULT_STENCIL_FRACTION);

    template<scalar::Field ScalarFieldType>
    [[nodiscard]] Result interpolate(const TriangleMesh& mesh, ScalarFieldType& scalar_field) const;

 private:
    // The domain points of the six quadratics, deduplicated: the seven
    // vertices of the split, then the midpoint of each of its twelve edges.
    static constexpr int DOMAIN_COUNT = 19;
    static constexpr int VERTEX_CONDITION_COUNT = 12;
    static constexpr int SMOOTHNESS_CONDITION_COUNT = 12;
    static constexpr int DEGREE = 2;
    static constexpr double DEFAULT_STENCIL_FRACTION = 0.75;
    static constexpr int MAXIMUM_DAMPING_ROUNDS = 60;
    static constexpr double DAMPING_FACTOR = 0.5;

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

    // The domain points touching each corner. Setting them from the corner's
    // gradient is what pins the value and the gradient there: for a quadratic
    // the coefficient beside a corner is half the gradient's reading of the
    // corner it leans towards.
    static constexpr std::array<std::array<int, 4>, 3> CORNER_DOMAIN_POINTS = {{
        {0, 7, 12, 13},
        {1, 8, 9, 14},
        {2, 10, 11, 15}
    }};

    double _stencil_fraction;

    [[nodiscard]] static double mean_edge_length(const TriangleMesh& mesh);

    [[nodiscard]] static std::array<Vector3, 3> corner_gradients(
        const PowellSabinRefinement::Cell& cell,
        const std::vector<double>& values,
        const std::vector<Vector3>& gradients,
        const std::vector<double>& damping
    );

    [[nodiscard]] static Eigen::VectorXd solve_cell(
        const PowellSabinRefinement::Cell& cell,
        const std::array<Vector3, 3>& corner_gradients
    );

    [[nodiscard]] static std::array<Vector3, DOMAIN_COUNT> domain_points(const PowellSabinRefinement::Cell& cell);
    [[nodiscard]] static int domain_index(int first, int second);
    [[nodiscard]] static int domain_index_for(const std::array<int, 3>& corner_order, const MultiIndex& index);
    [[nodiscard]] static bool contains(const std::array<int, 3>& values, int value);
};

inline PowellSabinInterpolant::PowellSabinInterpolant(double stencil_fraction) :
    _stencil_fraction(stencil_fraction) {
}

template<scalar::Field ScalarFieldType>
PowellSabinInterpolant::Result PowellSabinInterpolant::interpolate(
    const TriangleMesh& mesh,
    ScalarFieldType& scalar_field
) const {
    PowellSabinRefinement refinement(mesh);

    LocalQuadraticFit reading(_stencil_fraction * mean_edge_length(mesh));
    std::vector<double> values(mesh.vertices.size());
    std::vector<Vector3> gradients(mesh.vertices.size());

    for (size_t index = 0; index < mesh.vertices.size(); ++index) {
        LocalQuadraticFit::Reading local = reading.at(scalar_field, mesh.vertices[index]);
        values[index] = local.value;
        gradients[index] = local.tangential_gradient;
    }

    const std::vector<PowellSabinRefinement::Cell>& cells = refinement.cells();
    std::vector<double> damping(mesh.vertices.size(), 1.0);
    std::vector<Eigen::VectorXd> coefficients(cells.size());
    std::vector<bool> stale(cells.size(), true);

    // Interpolating a gradient the mesh is too coarse to carry makes the
    // piece overshoot, and a density that dips below zero is not a density.
    // Damping is per vertex rather than per piece, so every piece meeting at
    // a vertex still sees the same data and the field stays C1 while it is
    // brought back into range.
    for (int round = 0; round <= MAXIMUM_DAMPING_ROUNDS; ++round) {
        for (size_t index = 0; index < cells.size(); ++index) {
            if (stale[index]) {
                coefficients[index] = solve_cell(
                    cells[index],
                    corner_gradients(cells[index], values, gradients, damping)
                );
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
        *std::min_element(damping.begin(), damping.end())
    };
}

// A piece is a quadratic that is homogeneous in space, so its radial
// derivative is twice its value: the value and the tangential gradient
// together name one vector, and every corner condition is a dot product
// against it.
inline std::array<Vector3, 3> PowellSabinInterpolant::corner_gradients(
    const PowellSabinRefinement::Cell& cell,
    const std::vector<double>& values,
    const std::vector<Vector3>& gradients,
    const std::vector<double>& damping
) {
    std::array<Vector3, 3> result;

    for (int corner = 0; corner < 3; ++corner) {
        size_t vertex = cell.corner_vertices[corner];
        result[corner] = 2.0 * values[vertex] * cell.points[corner] + damping[vertex] * gradients[vertex];
    }

    return result;
}

inline double PowellSabinInterpolant::mean_edge_length(const TriangleMesh& mesh) {
    double total = 0.0;
    size_t count = 0;

    for (const std::array<size_t, 3>& triangle : mesh.triangles) {
        for (int corner = 0; corner < 3; ++corner) {
            const VectorS2& from = mesh.vertices[triangle[corner]];
            const VectorS2& to = mesh.vertices[triangle[(corner + 1) % 3]];
            total += std::acos(std::clamp(from.dot(to), -1.0, 1.0));
            ++count;
        }
    }

    return total / static_cast<double>(count);
}

// One small dense solve per cell. The unknowns are the coefficients at the
// nineteen domain points, shared between neighbouring pieces so that
// continuity of value holds by construction; the rows say what the corners
// interpolate and that the pieces meet smoothly inside the cell. Smoothness
// across a cell's boundary comes from the refinement's geometry, not from
// these rows.
inline Eigen::VectorXd PowellSabinInterpolant::solve_cell(
    const PowellSabinRefinement::Cell& cell,
    const std::array<Vector3, 3>& corner_gradients
) {
    constexpr int ROW_COUNT = VERTEX_CONDITION_COUNT + SMOOTHNESS_CONDITION_COUNT;

    std::array<Vector3, DOMAIN_COUNT> points = domain_points(cell);
    Eigen::MatrixXd system = Eigen::MatrixXd::Zero(ROW_COUNT, DOMAIN_COUNT);
    Eigen::VectorXd values = Eigen::VectorXd::Zero(ROW_COUNT);
    int row = 0;

    for (int corner = 0; corner < 3; ++corner) {
        for (int domain : CORNER_DOMAIN_POINTS[corner]) {
            system(row, domain) = 1.0;
            values(row) = corner_gradients[corner].dot(points[domain] - 0.5 * cell.points[corner]);
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

    return system.colPivHouseholderQr().solve(values);
}

inline std::array<Vector3, PowellSabinInterpolant::DOMAIN_COUNT>
PowellSabinInterpolant::domain_points(const PowellSabinRefinement::Cell& cell) {
    std::array<Vector3, DOMAIN_COUNT> result;

    for (int point = 0; point < PowellSabinRefinement::POINT_COUNT; ++point) {
        result[point] = cell.points[point];
    }

    for (const std::array<int, 3>& midpoint : EDGE_MIDPOINTS) {
        result[midpoint[2]] = 0.5 * (cell.points[midpoint[0]] + cell.points[midpoint[1]]);
    }

    return result;
}

inline int PowellSabinInterpolant::domain_index(int first, int second) {
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

inline int PowellSabinInterpolant::domain_index_for(const std::array<int, 3>& corner_order, const MultiIndex& index) {
    std::array<int, DEGREE> points{};
    int count = 0;

    for (int axis = 0; axis < 3; ++axis) {
        for (int repeat = 0; repeat < index[axis]; ++repeat) {
            points[count++] = corner_order[axis];
        }
    }

    return domain_index(points[0], points[1]);
}

inline bool PowellSabinInterpolant::contains(const std::array<int, 3>& values, int value) {
    return std::find(values.begin(), values.end(), value) != values.end();
}

} // namespace globe::fields::spherical

#endif //GLOBEART_SRC_GLOBE_FIELDS_SPHERICAL_POWELL_SABIN_INTERPOLANT_HPP_
