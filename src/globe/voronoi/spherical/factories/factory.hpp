#ifndef GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_FACTORIES_FACTORY_HPP_
#define GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_FACTORIES_FACTORY_HPP_

#include "../core/sphere.hpp"
#include "../core/random_builder.hpp"
#include "../core/callback.hpp"
#include "../optimizers/field_density_optimizer.hpp"
#include "../optimizers/gradient_density_optimizer.hpp"
#include "../optimizers/lloyd_optimizer.hpp"
#include "../../../fields/scalar/noise_field.hpp"
#include "../../../fields/spherical/polynomial_field.hpp"
#include "../../../fields/spherical/polynomial_field_fitter.hpp"
#include "../../../generators/spherical/fibonacci_point_generator.hpp"
#include "../../../generators/spherical/random_point_generator.hpp"
#include <string>
#include <memory>
#include <algorithm>
#include <iostream>
#include <iomanip>

namespace globe::voronoi::spherical {

using fields::scalar::NoiseField;
using fields::spherical::PolynomialField;
using fields::spherical::PolynomialFieldFitter;

class Factory {
 public:
    Factory(
        int points_count,
        std::string density_function,
        std::string optimization_strategy,
        int optimization_passes,
        int lloyd_passes,
        int max_perturbations = 50,
        Callback callback = noop_callback()
    );

    std::unique_ptr<Sphere> build();

 private:
    static constexpr int NOISE_FIT_DEGREE = 8;
    static constexpr size_t NOISE_FIT_SAMPLES = 20000;

    int _points_count;
    std::string _density_function;
    std::string _optimization_strategy;
    size_t _optimization_passes;
    size_t _lloyd_passes;
    size_t _max_perturbations;
    Callback _callback;

    std::unique_ptr<Sphere> build_initial();
    std::unique_ptr<Sphere> optimize_density(std::unique_ptr<Sphere> sphere);

    PolynomialField create_field() const;
    std::unique_ptr<Sphere> optimize_ccvd(std::unique_ptr<Sphere> sphere, const PolynomialField& field);
    std::unique_ptr<Sphere> optimize_gradient(std::unique_ptr<Sphere> sphere, const PolynomialField& field);

    static PolynomialField fit_noise_field();
};

inline Factory::Factory(
    int points_count,
    std::string density_function,
    std::string optimization_strategy,
    int optimization_passes,
    int lloyd_passes,
    int max_perturbations,
    Callback callback
) :
    _points_count(points_count),
    _density_function(std::move(density_function)),
    _optimization_strategy(std::move(optimization_strategy)),
    _optimization_passes(static_cast<size_t>(optimization_passes)),
    _lloyd_passes(static_cast<size_t>(lloyd_passes)),
    _max_perturbations(static_cast<size_t>(max_perturbations)),
    _callback(std::move(callback)) {
}

inline std::unique_ptr<Sphere> Factory::build() {
    std::cout << "Generating " << _points_count << " random points..." << std::flush;
    auto sphere = build_initial();
    std::cout << " done" << std::endl;

    _callback(*sphere);

    std::cout << "Initial density optimization..." << std::endl;
    sphere = optimize_density(std::move(sphere));

    for (size_t i = 0; i < _lloyd_passes; i++) {
        LloydOptimizer lloyd_optimizer(std::move(sphere), 1, _callback);
        sphere = lloyd_optimizer.optimize();
        double deviation = lloyd_optimizer.final_deviation();

        std::cout << "  " << std::setw(8) << std::left << "Lloyd" << std::right <<
            std::setw(4) << (i + 1) << "/" <<
            std::setw(4) << std::left << _lloyd_passes << std::right <<
            ": dev " << std::fixed << std::setprecision(8) << deviation <<
            std::defaultfloat << std::endl;

        sphere = optimize_density(std::move(sphere));
    }

    return sphere;
}

inline PolynomialField Factory::create_field() const {
    if (_density_function == "constant") {
        return PolynomialField::constant(1.0);
    }

    if (_density_function == "linear") {
        return PolynomialField::linear(2.0, Eigen::Vector3d(0.0, 0.0, 2.0));
    }

    if (_density_function == "noise") {
        return fit_noise_field();
    }

    Eigen::Matrix3d quadratic = Eigen::Matrix3d::Zero();
    quadratic(2, 2) = -0.9;
    return PolynomialField::quadratic(1.0, Eigen::Vector3d::Zero(), quadratic);
}

inline PolynomialField Factory::fit_noise_field() {
    NoiseField noise_field;
    PolynomialFieldFitter<> fitter(NOISE_FIT_DEGREE, NOISE_FIT_SAMPLES, generators::spherical::FibonacciPointGenerator());
    auto fit = fitter.fit(noise_field);

    std::cout << "Fitted noise to degree " << NOISE_FIT_DEGREE <<
        " polynomial, RMS residual " << fit.root_mean_square_residual << std::endl;

    return fit.field;
}

inline std::unique_ptr<Sphere> Factory::optimize_density(
    std::unique_ptr<Sphere> sphere
) {
    PolynomialField field = create_field();

    if (_optimization_strategy == "gradient") {
        return optimize_gradient(std::move(sphere), field);
    }

    return optimize_ccvd(std::move(sphere), field);
}

inline std::unique_ptr<Sphere> Factory::build_initial() {
    RandomBuilder<> builder;
    return builder.build(_points_count);
}

inline std::unique_ptr<Sphere> Factory::optimize_ccvd(
    std::unique_ptr<Sphere> sphere,
    const PolynomialField& field
) {
    FieldDensityOptimizer optimizer(
        std::move(sphere),
        field,
        _optimization_passes,
        generators::spherical::RandomPointGenerator<>(),
        _callback
    );

    return optimizer.optimize();
}

inline std::unique_ptr<Sphere> Factory::optimize_gradient(
    std::unique_ptr<Sphere> sphere,
    const PolynomialField& field
) {
    GradientDensityOptimizer<PolynomialField> optimizer(
        std::move(sphere),
        field,
        _optimization_passes,
        _max_perturbations,
        generators::spherical::RandomPointGenerator<>(),
        _callback
    );

    return optimizer.optimize();
}

} // namespace globe::voronoi::spherical

#endif //GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_FACTORIES_FACTORY_HPP_
