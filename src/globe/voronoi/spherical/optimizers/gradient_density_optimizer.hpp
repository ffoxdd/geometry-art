#ifndef GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_GRADIENT_DENSITY_OPTIMIZER_HPP_
#define GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_GRADIENT_DENSITY_OPTIMIZER_HPP_

#include "../../../types.hpp"
#include "../../../geometry/spherical/helpers.hpp"
#include "../../../geometry/spherical/arc.hpp"
#include "../../../geometry/spherical/polygon/polygon.hpp"
#include "../../../fields/spherical/field.hpp"
#include "../../../fields/spherical/polynomial_field.hpp"
#include "../../../generators/spherical/point_generator.hpp"
#include "../../../generators/spherical/random_point_generator.hpp"
#include "../core/sphere.hpp"
#include "../core/callback.hpp"
#include <Eigen/Core>
#include <LBFGS.h>
#include <memory>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <optional>

namespace globe::voronoi::spherical {

template<
    fields::spherical::Field FieldType = fields::spherical::PolynomialField,
    generators::spherical::PointGenerator GeneratorType = generators::spherical::RandomPointGenerator<>
>
class GradientDensityOptimizer {
 public:
    static constexpr size_t DEFAULT_ITERATIONS = 100;
    static constexpr size_t DEFAULT_MAX_PERTURBATIONS = 50;
    static constexpr double CONVERGENCE_THRESHOLD = 1e-8;
    static constexpr size_t LBFGS_HISTORY_SIZE = 6;

    static constexpr double PERTURBATION_ANGLE = 0.05;
    static constexpr size_t MAX_PERTURBATIONS_BEFORE_RESTORE = 5;
    static constexpr size_t MAX_RESTORES_PER_CHECKPOINT = 3;

    GradientDensityOptimizer(
        std::unique_ptr<Sphere> voronoi_sphere,
        FieldType field,
        size_t max_iterations = DEFAULT_ITERATIONS,
        size_t max_perturbations = DEFAULT_MAX_PERTURBATIONS,
        GeneratorType point_generator = GeneratorType(),
        Callback callback = noop_callback()
    );

    std::unique_ptr<Sphere> optimize();

 private:
    using Checkpoint = std::vector<cgal::Point3>;

    class ObjectiveFunctor {
     public:
        ObjectiveFunctor(GradientDensityOptimizer& optimizer, size_t start_iteration);
        double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad);

        size_t iteration_count() const { return _iteration_count; }
        double last_rms_error() const { return _last_rms_error; }

     private:
        static constexpr size_t LOG_INTERVAL = 10;

        GradientDensityOptimizer& _optimizer;
        size_t _start_iteration;
        size_t _iteration_count = 0;
        double _last_rms_error = 0.0;
    };

    std::unique_ptr<Sphere> _sphere;
    FieldType _field;
    size_t _max_iterations;
    size_t _max_perturbations;
    double _target_mass;
    GeneratorType _point_generator;
    Callback _callback;

    Eigen::VectorXd sites_to_vector() const;
    void vector_to_sites_normalized(const Eigen::VectorXd& x);

    std::vector<double> compute_mass_errors() const;
    std::vector<Eigen::Vector3d> compute_gradients(
        const std::vector<double>& mass_errors
    ) const;
    Eigen::VectorXd compute_projected_gradient() const;
    double compute_total_error(const std::vector<double>& mass_errors) const;
    double compute_error() const;

    int run_lbfgs_phase(size_t max_iterations, size_t start_iteration);

    Checkpoint save_checkpoint() const;
    void restore_checkpoint(const Checkpoint& checkpoint);
    std::optional<size_t> find_most_undersized_index() const;
    bool perturb_site_randomly(size_t index);
};

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
GradientDensityOptimizer<FieldType, GeneratorType>::GradientDensityOptimizer(
    std::unique_ptr<Sphere> voronoi_sphere,
    FieldType field,
    size_t max_iterations,
    size_t max_perturbations,
    GeneratorType point_generator,
    Callback callback
) :
    _sphere(std::move(voronoi_sphere)),
    _field(std::move(field)),
    _max_iterations(max_iterations),
    _max_perturbations(max_perturbations),
    _target_mass(0.0),
    _point_generator(std::move(point_generator)),
    _callback(std::move(callback)) {
}

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
GradientDensityOptimizer<FieldType, GeneratorType>::ObjectiveFunctor::ObjectiveFunctor(
    GradientDensityOptimizer& optimizer,
    size_t start_iteration
) : _optimizer(optimizer), _start_iteration(start_iteration) {}

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
double GradientDensityOptimizer<FieldType, GeneratorType>::ObjectiveFunctor::operator()(
    const Eigen::VectorXd& x,
    Eigen::VectorXd& grad
) {
    _optimizer.vector_to_sites_normalized(x);

    double error = _optimizer.compute_error();
    grad = _optimizer.compute_projected_gradient();

    _last_rms_error = std::sqrt(2.0 * error / _optimizer._sphere->size());
    ++_iteration_count;
    _optimizer._callback(*_optimizer._sphere);

    return error;
}

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
int GradientDensityOptimizer<FieldType, GeneratorType>::run_lbfgs_phase(
    size_t max_iterations,
    size_t start_iteration
) {
    LBFGSpp::LBFGSParam<double> param;
    param.m = LBFGS_HISTORY_SIZE;
    param.epsilon = CONVERGENCE_THRESHOLD;
    param.max_iterations = static_cast<int>(max_iterations);
    param.past = 3;
    param.delta = 1e-8;

    LBFGSpp::LBFGSSolver<double> solver(param);
    ObjectiveFunctor functor(*this, start_iteration);

    Eigen::VectorXd x = sites_to_vector();
    double fx;

    try {
        int iterations = solver.minimize(functor, x, fx);
        vector_to_sites_normalized(x);
        return iterations;
    } catch (const std::exception&) {
        vector_to_sites_normalized(x);
        return -1;
    }
} // namespace globe::voronoi::spherical

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
std::unique_ptr<Sphere> GradientDensityOptimizer<FieldType, GeneratorType>::optimize() {
    _target_mass = _field.total_mass() / _sphere->size();

    double best_error = std::numeric_limits<double>::max();
    Checkpoint best_checkpoint = save_checkpoint();
    size_t perturbation_attempts = 0;
    size_t perturbations_since_best = 0;
    size_t restores_to_current_best = 0;
    size_t total_iterations = 0;
    size_t stalled_phases = 0;
    static constexpr size_t MAX_STALLED_PHASES = 3;

    while (total_iterations < _max_iterations) {
        size_t remaining = _max_iterations - total_iterations;
        size_t phase_iterations = std::min(remaining, static_cast<size_t>(100));
        size_t start_iteration = total_iterations;

        int result = run_lbfgs_phase(phase_iterations, start_iteration);
        total_iterations += (result > 0) ? static_cast<size_t>(result) : phase_iterations;

        _callback(*_sphere);

        double error = compute_error();
        double rms_error = std::sqrt(2.0 * error / _sphere->size());

        std::cout << "  " << std::setw(8) << std::left << "L-BFGS" << std::right <<
            std::setw(4) << start_iteration << "-" <<
            std::setw(4) << std::left << total_iterations << std::right <<
            ": RMS " << std::fixed << std::setprecision(8) << rms_error << std::defaultfloat;

        if (rms_error < CONVERGENCE_THRESHOLD) {
            std::cout << " [converged]" << std::endl;
            break;
        }

        bool is_new_best = error < best_error;
        bool significant_improvement = error < best_error * 0.999;

        if (is_new_best) {
            best_error = error;
            best_checkpoint = save_checkpoint();
            perturbations_since_best = 0;
            restores_to_current_best = 0;
            stalled_phases = 0;
        } else {
            ++stalled_phases;
        }

        bool should_perturb = !significant_improvement &&
                              (rms_error > CONVERGENCE_THRESHOLD * 100) &&
                              (result != 0);

        if (!should_perturb && stalled_phases >= MAX_STALLED_PHASES) {
            std::cout << " [stalled]" << std::endl;
            break;
        }

        if (should_perturb) {
            ++perturbations_since_best;

            if (perturbations_since_best >= MAX_PERTURBATIONS_BEFORE_RESTORE &&
                restores_to_current_best < MAX_RESTORES_PER_CHECKPOINT) {
                restore_checkpoint(best_checkpoint);
                ++restores_to_current_best;
                perturbations_since_best = 0;
                std::cout << " [restore #" << restores_to_current_best << "]" << std::endl;
                continue;
            }

            if (perturbation_attempts >= _max_perturbations) {
                std::cout << " [max perturb]" << std::endl;
                break;
            }

            auto undersized = find_most_undersized_index();
            if (undersized.has_value() && perturb_site_randomly(undersized.value())) {
                ++perturbation_attempts;
                if (is_new_best) {
                    std::cout << " [new best, perturb #" << perturbation_attempts << "]" << std::endl;
                } else {
                    std::cout << " [perturb #" << perturbation_attempts << "]" << std::endl;
                }
            } else {
                if (is_new_best) {
                    std::cout << " [new best]" << std::endl;
                } else {
                    std::cout << std::endl;
                }
            }
        } else if (is_new_best) {
            std::cout << " [new best]" << std::endl;
        } else {
            std::cout << std::endl;
        }
    }

    restore_checkpoint(best_checkpoint);

    double final_error = compute_error();
    double final_rms = std::sqrt(2.0 * final_error / _sphere->size());
    std::cout << "  " << std::setw(8) << std::left << "Final" << std::right <<
        std::setw(9) << "" <<
        ": RMS " << std::fixed << std::setprecision(8) << final_rms << std::defaultfloat << std::endl;

    return std::move(_sphere);
} // namespace globe::voronoi::spherical

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
Eigen::VectorXd GradientDensityOptimizer<FieldType, GeneratorType>::sites_to_vector() const {
    size_t n = _sphere->size();
    Eigen::VectorXd x(3 * n);

    for (size_t i = 0; i < n; ++i) {
        cgal::Point3 site = _sphere->site(i);
        x[3 * i + 0] = site.x();
        x[3 * i + 1] = site.y();
        x[3 * i + 2] = site.z();
    }

    return x;
} // namespace globe::voronoi::spherical

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
void GradientDensityOptimizer<FieldType, GeneratorType>::vector_to_sites_normalized(const Eigen::VectorXd& x) {
    size_t n = _sphere->size();

    for (size_t i = 0; i < n; ++i) {
        Eigen::Vector3d v(x[3 * i + 0], x[3 * i + 1], x[3 * i + 2]);
        v.normalize();
        _sphere->update_site(i, cgal::Point3(v.x(), v.y(), v.z()));
    }
} // namespace globe::voronoi::spherical

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
double GradientDensityOptimizer<FieldType, GeneratorType>::compute_error() const {
    auto mass_errors = compute_mass_errors();
    return compute_total_error(mass_errors);
} // namespace globe::voronoi::spherical

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
Eigen::VectorXd GradientDensityOptimizer<FieldType, GeneratorType>::compute_projected_gradient() const {
    auto mass_errors = compute_mass_errors();
    auto gradients = compute_gradients(mass_errors);

    size_t n = _sphere->size();
    Eigen::VectorXd grad(3 * n);

    for (size_t i = 0; i < n; ++i) {
        cgal::Point3 site = _sphere->site(i);
        Eigen::Vector3d s(site.x(), site.y(), site.z());
        s.normalize();

        Eigen::Vector3d g = gradients[i];
        g -= g.dot(s) * s;

        grad[3 * i + 0] = g.x();
        grad[3 * i + 1] = g.y();
        grad[3 * i + 2] = g.z();
    }

    return grad;
} // namespace globe::voronoi::spherical

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
std::vector<double> GradientDensityOptimizer<FieldType, GeneratorType>::compute_mass_errors() const {
    size_t n = _sphere->size();
    std::vector<double> errors(n);

    size_t i = 0;
    for (const auto& cell : _sphere->cells()) {
        errors[i] = _field.integrals(cell).mass - _target_mass;
        ++i;
    }

    return errors;
} // namespace globe::voronoi::spherical

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
std::vector<Eigen::Vector3d> GradientDensityOptimizer<FieldType, GeneratorType>::compute_gradients(
    const std::vector<double>& mass_errors
) const {
    size_t n = _sphere->size();
    std::vector<Eigen::Vector3d> gradients(n, Eigen::Vector3d::Zero());

    for (size_t k = 0; k < n; ++k) {
        cgal::Point3 site_k = _sphere->site(k);
        Eigen::Vector3d s_k = to_eigen(site_k);
        auto cell_edges = _sphere->cell_edges(k);

        for (const auto& edge_info : cell_edges) {
            size_t j = edge_info.neighbor_index;
            cgal::Point3 site_j = _sphere->site(j);

            auto edge = _field.integrals(edge_info.arc);
            Eigen::Vector3d rho_weighted_moment = edge.first_moment;
            double edge_integral = edge.mass;

            Eigen::Vector3d n_vec = to_eigen(site_j) - s_k;
            double n_norm = n_vec.norm();

            if (n_norm < GEOMETRIC_EPSILON) {
                continue;
            }

            Eigen::Vector3d edge_grad = (rho_weighted_moment - s_k * edge_integral) / n_norm;
            gradients[k] += (mass_errors[k] - mass_errors[j]) * edge_grad;
        }
    }

    return gradients;
} // namespace globe::voronoi::spherical

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
double GradientDensityOptimizer<FieldType, GeneratorType>::compute_total_error(
    const std::vector<double>& mass_errors
) const {
    double sum = 0.0;
    for (double e : mass_errors) {
        sum += e * e;
    }
    return sum / 2.0;
} // namespace globe::voronoi::spherical

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
typename GradientDensityOptimizer<FieldType, GeneratorType>::Checkpoint
GradientDensityOptimizer<FieldType, GeneratorType>::save_checkpoint() const {
    Checkpoint checkpoint;
    checkpoint.reserve(_sphere->size());

    for (size_t i = 0; i < _sphere->size(); ++i) {
        checkpoint.push_back(_sphere->site(i));
    }

    return checkpoint;
} // namespace globe::voronoi::spherical

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
void GradientDensityOptimizer<FieldType, GeneratorType>::restore_checkpoint(const Checkpoint& checkpoint) {
    auto new_voronoi = std::make_unique<Sphere>();
    for (const auto& site : checkpoint) {
        new_voronoi->insert(site);
    }
    _sphere = std::move(new_voronoi);
} // namespace globe::voronoi::spherical

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
std::optional<size_t> GradientDensityOptimizer<FieldType, GeneratorType>::find_most_undersized_index() const {
    auto mass_errors = compute_mass_errors();

    size_t best_index = 0;
    double largest_deficit = 0.0;
    double largest_abs_error = 0.0;
    size_t largest_error_index = 0;

    for (size_t i = 0; i < mass_errors.size(); ++i) {
        double deficit = -mass_errors[i];
        if (deficit > largest_deficit) {
            largest_deficit = deficit;
            best_index = i;
        }
        double abs_error = std::abs(mass_errors[i]);
        if (abs_error > largest_abs_error) {
            largest_abs_error = abs_error;
            largest_error_index = i;
        }
    }

    if (largest_deficit > 0.0) {
        return best_index;
    }

    if (largest_abs_error > 0.0) {
        return largest_error_index;
    }

    return std::nullopt;
} // namespace globe::voronoi::spherical

template<fields::spherical::Field FieldType, generators::spherical::PointGenerator GeneratorType>
bool GradientDensityOptimizer<FieldType, GeneratorType>::perturb_site_randomly(size_t index) {
    cgal::Point3 site = _sphere->site(index);
    Eigen::Vector3d s(site.x(), site.y(), site.z());
    s.normalize();

    auto random_points = _point_generator.generate(1);
    if (random_points.empty()) {
        return false;
    }

    VectorS2 random_point = random_points[0];
    Eigen::Vector3d random_dir = random_point;

    Eigen::Vector3d tangent = random_dir - random_dir.dot(s) * s;
    double tangent_norm = tangent.norm();

    if (tangent_norm < GEOMETRIC_EPSILON) {
        return false;
    }

    tangent /= tangent_norm;

    Eigen::Vector3d new_pos =
        std::cos(PERTURBATION_ANGLE) * s +
        std::sin(PERTURBATION_ANGLE) * tangent;

    new_pos.normalize();
    _sphere->update_site(index, cgal::Point3(new_pos.x(), new_pos.y(), new_pos.z()));

    return true;
} // namespace globe::voronoi::spherical

} // namespace globe::voronoi::spherical

#endif //GLOBEART_SRC_GLOBE_VORONOI_SPHERICAL_OPTIMIZERS_GRADIENT_DENSITY_OPTIMIZER_HPP_
