#ifndef GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_FACTORIES_FACTORY_HPP_
#define GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_FACTORIES_FACTORY_HPP_

#include "../core/sphere.hpp"
#include "../core/random_builder.hpp"
#include "../core/callback.hpp"
#include "../optimizers/capacity_constrained_optimizer.hpp"
#include "../optimizers/lloyd_optimizer.hpp"
#include "../../../fields/scalar/noise_field.hpp"
#include "../../../fields/spherical/polynomial_field.hpp"
#include "../../../fields/spherical/polynomial_field_fitter.hpp"
#include "../../../generators/spherical/fibonacci_point_generator.hpp"
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <utility>

namespace globe::voronoi::spherical {

using fields::scalar::NoiseField;
using fields::spherical::PolynomialField;
using fields::spherical::PolynomialFieldFitter;

class Factory {
 public:
    Factory(
        int points_count,
        std::string density_field,
        size_t lloyd_passes,
        CapacityConstrainedParameters optimizer_parameters,
        Callback callback
    );

    std::unique_ptr<Sphere> build();

 private:
    static constexpr int NOISE_FIT_DEGREE = 8;
    static constexpr size_t NOISE_FIT_SAMPLES = 20000;

    int _points_count;
    std::string _density_field;
    size_t _lloyd_passes;
    CapacityConstrainedParameters _optimizer_parameters;
    Callback _callback;

    [[nodiscard]] PolynomialField create_field() const;
    [[nodiscard]] std::unique_ptr<Sphere> build_initial() const;
    [[nodiscard]] std::unique_ptr<Sphere> warm_start(std::unique_ptr<Sphere> sphere, const PolynomialField& field) const;
    [[nodiscard]] std::unique_ptr<Sphere> optimize(std::unique_ptr<Sphere> sphere, const PolynomialField& field) const;

    [[nodiscard]] static PolynomialField fit_noise_field();
};

inline Factory::Factory(
    int points_count,
    std::string density_field,
    size_t lloyd_passes,
    CapacityConstrainedParameters optimizer_parameters,
    Callback callback
) :
    _points_count(points_count),
    _density_field(std::move(density_field)),
    _lloyd_passes(lloyd_passes),
    _optimizer_parameters(optimizer_parameters),
    _callback(std::move(callback)) {
}

inline std::unique_ptr<Sphere> Factory::build() {
    PolynomialField field = create_field();

    std::cout << "Generating " << _points_count << " random points..." << std::flush;
    auto sphere = build_initial();
    std::cout << " done" << std::endl;
    _callback(*sphere);

    sphere = warm_start(std::move(sphere), field);
    return optimize(std::move(sphere), field);
}

inline PolynomialField Factory::create_field() const {
    if (_density_field == "constant") {
        return PolynomialField::constant(1.0);
    }

    if (_density_field == "linear") {
        return PolynomialField::linear(2.0, Vector3(0.0, 0.0, 2.0));
    }

    if (_density_field == "noise") {
        return fit_noise_field();
    }

    Eigen::Matrix3d quadratic = Eigen::Matrix3d::Zero();
    quadratic(2, 2) = -0.9;
    return PolynomialField::quadratic(1.0, Vector3::Zero(), quadratic);
}

inline PolynomialField Factory::fit_noise_field() {
    NoiseField noise_field;
    PolynomialFieldFitter<> fitter(NOISE_FIT_DEGREE, NOISE_FIT_SAMPLES, generators::spherical::FibonacciPointGenerator());
    auto fit = fitter.fit(noise_field);

    std::cout << "Fitted noise to degree " << NOISE_FIT_DEGREE <<
        " polynomial, RMS residual " << fit.root_mean_square_residual << std::endl;

    return fit.field;
}

inline std::unique_ptr<Sphere> Factory::build_initial() const {
    return RandomBuilder<>().build(_points_count);
}

inline std::unique_ptr<Sphere> Factory::warm_start(std::unique_ptr<Sphere> sphere, const PolynomialField& field) const {
    if (_lloyd_passes == 0) {
        return sphere;
    }

    LloydOptimizer<PolynomialField> lloyd(std::move(sphere), field, _lloyd_passes, _callback);
    sphere = lloyd.optimize();

    std::cout << "  " << std::setw(8) << std::left << "Lloyd" << std::right <<
        std::setw(3) << _lloyd_passes << " passes: centroid deviation " <<
        std::scientific << std::setprecision(3) << lloyd.final_deviation() <<
        std::defaultfloat << std::endl;

    return sphere;
}

inline std::unique_ptr<Sphere> Factory::optimize(std::unique_ptr<Sphere> sphere, const PolynomialField& field) const {
    CapacityConstrainedOptimizer<PolynomialField> optimizer(std::move(sphere), field, _optimizer_parameters, _callback);
    sphere = optimizer.optimize();

    const auto& report = optimizer.report();
    std::cout << "  " << std::setw(8) << std::left << "Final" << std::right <<
        (report.converged ? " converged" : (report.stalled ? " stalled" : " stopped")) <<
        " after " << report.outer_iterations << " outer / " << report.inner_iterations << " inner iterations" <<
        ", capacity RMS " << std::scientific << std::setprecision(3) << report.relative_rms_capacity_error <<
        std::defaultfloat << std::endl;

    return sphere;
}

} // namespace globe::voronoi::spherical

#endif //GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_FACTORIES_FACTORY_HPP_
