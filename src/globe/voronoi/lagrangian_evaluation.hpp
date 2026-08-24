#ifndef GLOBEART_SRC_GLOBE_VORONOI_LAGRANGIAN_EVALUATION_HPP_
#define GLOBEART_SRC_GLOBE_VORONOI_LAGRANGIAN_EVALUATION_HPP_

#include "../types.hpp"
#include <cmath>
#include <cstddef>
#include <vector>

namespace globe::voronoi {

struct LagrangianEvaluation {
    double value;
    double cvt_energy;
    std::vector<double> capacity_errors;
    std::vector<Vector3> site_gradients;

    [[nodiscard]] double root_mean_square_capacity_error() const;
    [[nodiscard]] double max_absolute_capacity_error() const;
};

inline double LagrangianEvaluation::root_mean_square_capacity_error() const {
    double sum = 0.0;

    for (double error : capacity_errors) {
        sum += error * error;
    }

    return std::sqrt(sum / static_cast<double>(capacity_errors.size()));
}

inline double LagrangianEvaluation::max_absolute_capacity_error() const {
    double maximum = 0.0;

    for (double error : capacity_errors) {
        maximum = std::max(maximum, std::abs(error));
    }

    return maximum;
}

} // namespace globe::voronoi

#endif //GLOBEART_SRC_GLOBE_VORONOI_LAGRANGIAN_EVALUATION_HPP_
