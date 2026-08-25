#ifndef GEOMETRY_ART_FIELDS_SCALAR_NOISE_FIELD_HPP_
#define GEOMETRY_ART_FIELDS_SCALAR_NOISE_FIELD_HPP_

#include "../../types.hpp"
#include "../../math/interval.hpp"
#include <anl/anl.h>
#include <cmath>

namespace geometry_art::fields::scalar {

class NoiseField {
 public:
    NoiseField(Interval output_range = Interval(0, 1), int seed = 1546);
    double value(const VectorS2 &point);
    double max_frequency() const;

 private:
    anl::CKernel _kernel;
    anl::CInstructionIndex _instruction_index;

    static constexpr double NOISE_FREQUENCY = 1.0;
    static constexpr int NOISE_OCTAVES = 2;
    static constexpr double NOISE_PERSISTENCE = 0.5;
    static constexpr double NOISE_LACUNARITY = 2.0;
    static constexpr double CONTRAST_GAIN = 0.95;

    static double amplitude_sum(double persistence, int octaves);
    static double fractal_range_scale(double persistence, int octaves);
    static anl::CInstructionIndex initialize_kernel(anl::CKernel &kernel, Interval output_range, int seed);
};

inline NoiseField::NoiseField(Interval output_range, int seed) :
    _kernel(),
    _instruction_index(initialize_kernel(_kernel, output_range, seed)) {
}

inline double NoiseField::value(const VectorS2 &point) {
    return anl::CNoiseExecutor(_kernel).evaluateScalar(
        point.x(), point.y(), point.z(), _instruction_index
    );
}

inline double NoiseField::max_frequency() const {
    return NOISE_FREQUENCY * std::pow(NOISE_LACUNARITY, NOISE_OCTAVES - 1);
}

inline double NoiseField::amplitude_sum(double persistence, int octaves) {
    double sum = 0.0;

    for (int i = 0; i < octaves; ++i) {
        sum += std::pow(persistence, i);
    }

    return sum;
}

inline double NoiseField::fractal_range_scale(double persistence, int octaves) {
    constexpr double SIMPLEX_EMPIRICAL_FACTOR = 10.5;

    double amplitude = amplitude_sum(persistence, octaves);
    double range = 2.0 * SIMPLEX_EMPIRICAL_FACTOR * amplitude;

    return 1.0 / range;
}

inline anl::CInstructionIndex NoiseField::initialize_kernel(anl::CKernel &kernel, Interval output_range, int seed) {
    auto seed_instruction = kernel.constant(seed);

    auto instruction_index = kernel.fractal(
        seed_instruction,
        kernel.simplexBasis(seed_instruction),
        kernel.constant(NOISE_PERSISTENCE),
        kernel.constant(NOISE_LACUNARITY),
        kernel.constant(NOISE_OCTAVES),
        kernel.constant(NOISE_FREQUENCY)
    );

    double scale = fractal_range_scale(NOISE_PERSISTENCE, NOISE_OCTAVES);
    instruction_index = kernel.scaleOffset(instruction_index, scale, 0.5);
    instruction_index = kernel.gain(kernel.constant(CONTRAST_GAIN), instruction_index);

    double range_scale = output_range.measure();
    double range_offset = output_range.low();
    instruction_index = kernel.scaleOffset(instruction_index, range_scale, range_offset);

    instruction_index = kernel.clamp(
        instruction_index,
        kernel.constant(output_range.low()),
        kernel.constant(output_range.high())
    );

    return instruction_index;
}

} // namespace geometry_art::fields::scalar

#endif //GEOMETRY_ART_FIELDS_SCALAR_NOISE_FIELD_HPP_
