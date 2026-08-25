#ifndef GEOMETRY_ART_TESTING_ASSERTIONS_STATISTICAL_HPP_
#define GEOMETRY_ART_TESTING_ASSERTIONS_STATISTICAL_HPP_

#include <limits>
#include <vector>
#include <gtest/gtest.h>
#include "../../cgal/types.hpp"
#include "../../math/interval.hpp"

namespace geometry_art::testing::assertions {

struct DistributionMetrics {
    double mean;
    double variance;
    double min_value;
    double max_value;
    size_t sample_count;
    size_t values_at_min;
    size_t values_at_max;
    double range_span;
    double clipping_ratio;
};

struct CoordinateMetrics {
    DistributionMetrics x;
    DistributionMetrics y;
    DistributionMetrics z;
};

inline DistributionMetrics compute_statistics(const std::vector<double>& values) {
    size_t sample_count = values.size();
    double mean = 0.0;
    double m2 = 0.0;
    double min_value = std::numeric_limits<double>::max();
    double max_value = std::numeric_limits<double>::lowest();
    size_t values_at_min = 0;
    size_t values_at_max = 0;

    for (size_t i = 0; i < sample_count; ++i) {
        double value = values[i];

        double delta = value - mean;
        mean += delta / (i + 1);
        double delta2 = value - mean;
        m2 += delta * delta2;

        if (value < min_value) {
            min_value = value;
            values_at_min = 1;
        } else if (value == min_value) {
            values_at_min++;
        }

        if (value > max_value) {
            max_value = value;
            values_at_max = 1;
        } else if (value == max_value) {
            values_at_max++;
        }
    }

    double variance = m2 / sample_count;
    double range_span = max_value - min_value;
    double clipping_ratio = static_cast<double>(values_at_min + values_at_max) / sample_count;

    return {
        mean,
        variance,
        min_value,
        max_value,
        sample_count,
        values_at_min,
        values_at_max,
        range_span,
        clipping_ratio
    };
}

template<typename SampleFunc>
DistributionMetrics compute_statistics(
    SampleFunc sample_func,
    size_t sample_count
) {
    std::vector<double> values;
    values.reserve(sample_count);

    for (size_t i = 0; i < sample_count; ++i) {
        values.push_back(sample_func());
    }

    return compute_statistics(values);
}

inline CoordinateMetrics compute_coordinate_statistics(const std::vector<cgal::Point3>& points) {
    size_t sample_count = points.size();

    auto compute_from_coordinate = [&points, sample_count](auto coordinate_extractor) {
        size_t index = 0;
        return compute_statistics(
            [&]() { return coordinate_extractor(points[index++]); },
            sample_count
        );
    };

    return CoordinateMetrics{
        compute_from_coordinate([](const cgal::Point3& p) { return p.x(); }),
        compute_from_coordinate([](const cgal::Point3& p) { return p.y(); }),
        compute_from_coordinate([](const cgal::Point3& p) { return p.z(); })
    };
}

inline CoordinateMetrics compute_coordinate_statistics(const std::vector<VectorS2>& points) {
    size_t sample_count = points.size();

    auto compute_from_coordinate = [&points, sample_count](auto coordinate_extractor) {
        size_t index = 0;
        return compute_statistics(
            [&]() { return coordinate_extractor(points[index++]); },
            sample_count
        );
    };

    return CoordinateMetrics{
        compute_from_coordinate([](const VectorS2& p) { return p.x(); }),
        compute_from_coordinate([](const VectorS2& p) { return p.y(); }),
        compute_from_coordinate([](const VectorS2& p) { return p.z(); })
    };
}

template<typename PointFunc>
CoordinateMetrics compute_coordinate_statistics(
    PointFunc point_func,
    size_t sample_count
) {
    std::vector<cgal::Point3> points;
    points.reserve(sample_count);

    for (size_t i = 0; i < sample_count; ++i) {
        points.push_back(point_func());
    }

    return compute_coordinate_statistics(points);
}

inline void expect_mean(
    const DistributionMetrics& metrics,
    double expected_mean,
    double tolerance
) {
    EXPECT_NEAR(metrics.mean, expected_mean, tolerance)
        << "Expected mean " << expected_mean << " but got " << metrics.mean
        << " (sample count: " << metrics.sample_count << ")";
}

inline void expect_variance(
    const DistributionMetrics& metrics,
    double expected_variance,
    double tolerance
) {
    EXPECT_NEAR(metrics.variance, expected_variance, tolerance)
        << "Expected variance " << expected_variance << " but got " << metrics.variance
        << " (sample count: " << metrics.sample_count << ")";
}

inline double uniform_distribution_variance(double range_low, double range_high) {
    double range = range_high - range_low;
    return (range * range) / 12.0;
}

inline double uniform_distribution_mean(double range_low, double range_high) {
    return (range_low + range_high) / 2.0;
}

inline void expect_range_coverage(
    const DistributionMetrics& metrics,
    const Interval& expected_range,
    double expected_coverage,
    double tolerance = 0.0
) {
    double expected_span = expected_range.measure();
    double actual_coverage = 0.0;

    if (expected_span != 0.0) {
        actual_coverage = metrics.range_span / expected_span;
    }

    if (tolerance > 0.0) {
        EXPECT_NEAR(actual_coverage, expected_coverage, tolerance)
            << "Expected range coverage " << expected_coverage << " but got " << actual_coverage
            << " (range span: " << metrics.range_span << ", expected span: " << expected_span << ")";
    } else {
        EXPECT_GE(actual_coverage, expected_coverage)
            << "Expected range coverage >= " << expected_coverage << " but got " << actual_coverage
            << " (range span: " << metrics.range_span << ", expected span: " << expected_span << ")";
    }
}

} // namespace geometry_art::testing::assertions

namespace geometry_art::testing {
using assertions::DistributionMetrics;
using assertions::CoordinateMetrics;
using assertions::compute_statistics;
using assertions::compute_coordinate_statistics;
using assertions::expect_mean;
using assertions::expect_variance;
using assertions::uniform_distribution_variance;
using assertions::uniform_distribution_mean;
using assertions::expect_range_coverage;
}

#endif //GEOMETRY_ART_TESTING_ASSERTIONS_STATISTICAL_HPP_
