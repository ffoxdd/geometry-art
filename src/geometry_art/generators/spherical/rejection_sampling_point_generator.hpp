#ifndef GEOMETRY_ART_GENERATORS_SPHERICAL_REJECTION_SAMPLING_POINT_GENERATOR_HPP_
#define GEOMETRY_ART_GENERATORS_SPHERICAL_REJECTION_SAMPLING_POINT_GENERATOR_HPP_

#include "../../types.hpp"
#include "../../geometry/spherical/bounding_box.hpp"
#include "../../fields/scalar/field.hpp"
#include "../../fields/scalar/constant_field.hpp"
#include "../../math/interval_sampler/interval_sampler.hpp"
#include "point_generator.hpp"
#include "random_point_generator.hpp"
#include <vector>
#include <cmath>

namespace geometry_art::generators::spherical {

template<
    fields::scalar::Field DensityFieldType = fields::scalar::ConstantField,
    PointGenerator UnderlyingGeneratorType = RandomPointGenerator<>,
    IntervalSampler IntervalSamplerType = UniformIntervalSampler
>
class RejectionSamplingPointGenerator {
 public:
    RejectionSamplingPointGenerator() = default;

    explicit RejectionSamplingPointGenerator(
        DensityFieldType density_field,
        double max_density = 1.0,
        UnderlyingGeneratorType generator = UnderlyingGeneratorType(),
        IntervalSamplerType interval_sampler = IntervalSamplerType()
    ) : _density_field(std::move(density_field)),
        _max_density(std::max(max_density, 1e-6)),
        _generator(std::move(generator)),
        _interval_sampler(std::move(interval_sampler)) {
    }

    std::vector<VectorS2> generate(size_t count);
    std::vector<VectorS2> generate(
        size_t count,
        const SphericalBoundingBox &bounding_box
    );

    [[nodiscard]] size_t last_attempt_count() const { return _last_attempt_count; }

 private:
    DensityFieldType _density_field;
    double _max_density = 1.0;
    UnderlyingGeneratorType _generator;
    IntervalSamplerType _interval_sampler;
    size_t _last_attempt_count = 0;

    static constexpr size_t BATCH_SIZE = 100;

    [[nodiscard]] bool candidate_accepted(const VectorS2 &candidate);
    [[nodiscard]] double acceptance_threshold(const VectorS2 &point);
    [[nodiscard]] double sample_acceptance();
};

template<fields::scalar::Field DensityFieldType, PointGenerator UnderlyingGeneratorType, IntervalSampler IntervalSamplerType>
std::vector<VectorS2> RejectionSamplingPointGenerator<DensityFieldType, UnderlyingGeneratorType, IntervalSamplerType>::generate(
    size_t count
) {
    return generate(count, SphericalBoundingBox::full_sphere());
}

template<fields::scalar::Field DensityFieldType, PointGenerator UnderlyingGeneratorType, IntervalSampler IntervalSamplerType>
std::vector<VectorS2> RejectionSamplingPointGenerator<DensityFieldType, UnderlyingGeneratorType, IntervalSamplerType>::generate(
    size_t count,
    const SphericalBoundingBox &bounding_box
) {
    if (count == 0) {
        _last_attempt_count = 0;
        return {};
    }

    std::vector<VectorS2> points;
    points.reserve(count);
    size_t attempts = 0;

    while (points.size() < count) {
        size_t remaining = count - points.size();
        size_t batch_size = std::min(remaining, BATCH_SIZE);
        auto candidates = _generator.generate(batch_size, bounding_box);
        attempts += candidates.size();

        for (const auto &candidate : candidates) {
            if (points.size() >= count) {
                break;
            }

            if (candidate_accepted(candidate)) {
                points.push_back(candidate);
            }
        }
    }

    _last_attempt_count = attempts;
    return points;
}

template<fields::scalar::Field DensityFieldType, PointGenerator UnderlyingGeneratorType, IntervalSampler IntervalSamplerType>
bool RejectionSamplingPointGenerator<DensityFieldType, UnderlyingGeneratorType, IntervalSamplerType>::candidate_accepted(
    const VectorS2 &candidate
) {
    return sample_acceptance() <= acceptance_threshold(candidate);
}

template<fields::scalar::Field DensityFieldType, PointGenerator UnderlyingGeneratorType, IntervalSampler IntervalSamplerType>
double RejectionSamplingPointGenerator<DensityFieldType, UnderlyingGeneratorType, IntervalSamplerType>::acceptance_threshold(
    const VectorS2 &point
) {
    double density = _density_field.value(point);
    return std::clamp(density / _max_density, 0.0, 1.0);
}

template<fields::scalar::Field DensityFieldType, PointGenerator UnderlyingGeneratorType, IntervalSampler IntervalSamplerType>
double RejectionSamplingPointGenerator<DensityFieldType, UnderlyingGeneratorType, IntervalSamplerType>::sample_acceptance() {
    return _interval_sampler.sample(Interval(0.0, 1.0));
}

} // namespace geometry_art::generators::spherical

#endif //GEOMETRY_ART_GENERATORS_SPHERICAL_REJECTION_SAMPLING_POINT_GENERATOR_HPP_
