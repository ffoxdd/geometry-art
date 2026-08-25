#ifndef GEOMETRY_ART_FIELDS_SCALAR_IMAGE_FIELD_HPP_
#define GEOMETRY_ART_FIELDS_SCALAR_IMAGE_FIELD_HPP_

#include "../../types.hpp"
#include "../../math/interval.hpp"
#include <anl/Imaging/stb_image.h>
#include <CGAL/assertions.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace geometry_art::fields::scalar {

using geometry_art::math::Interval;

// A density read from an image in the equirectangular projection: longitude
// runs across, latitude runs down from the north pole, and darker pixels
// are denser, so the cells crowd where the ink is.
//
// The image is sampled bilinearly, wrapping in longitude and clamping at
// the poles. Structure finer than a cell cannot survive into the
// tessellation anyway, so no further smoothing is needed here: the spline
// projection averages it away.
class ImageField {
 public:
    ImageField(std::vector<double> brightness, size_t width, size_t height, Interval density_range);

    [[nodiscard]] static ImageField load(const std::string& path, Interval density_range);

    [[nodiscard]] double value(const VectorS2& point) const;

    // The density at a fraction of the picture: u runs across and wraps, v
    // runs down from the top and clamps.
    [[nodiscard]] double value_at(double u, double v) const;

    [[nodiscard]] size_t width() const { return _width; }
    [[nodiscard]] size_t height() const { return _height; }

 private:
    std::vector<double> _brightness;
    size_t _width;
    size_t _height;
    Interval _density_range;

    [[nodiscard]] double brightness_at(long row, long column) const;
};

inline ImageField::ImageField(
    std::vector<double> brightness,
    size_t width,
    size_t height,
    Interval density_range
) :
    _brightness(std::move(brightness)),
    _width(width),
    _height(height),
    _density_range(density_range) {
    CGAL_precondition(_brightness.size() == width * height);
    CGAL_precondition(width > 0 && height > 0);
}

inline ImageField ImageField::load(const std::string& path, Interval density_range) {
    int width = 0;
    int height = 0;
    int channels = 0;
    unsigned char* pixels = stbi_load(path.c_str(), &width, &height, &channels, 0);

    if (pixels == nullptr) {
        throw std::runtime_error("could not load image: " + path);
    }

    std::vector<double> brightness;
    brightness.reserve(static_cast<size_t>(width) * static_cast<size_t>(height));

    for (long pixel = 0; pixel < static_cast<long>(width) * height; ++pixel) {
        const unsigned char* sample = pixels + pixel * channels;
        double level = channels >= 3
            ? 0.2126 * sample[0] + 0.7152 * sample[1] + 0.0722 * sample[2]
            : static_cast<double>(sample[0]);
        brightness.push_back(level / 255.0);
    }

    stbi_image_free(pixels);

    return ImageField(
        std::move(brightness),
        static_cast<size_t>(width),
        static_cast<size_t>(height),
        density_range
    );
}

inline double ImageField::value(const VectorS2& point) const {
    double longitude = std::atan2(point.y(), point.x());
    double latitude = std::asin(std::clamp(point.z(), -1.0, 1.0));

    return value_at((longitude + M_PI) / (2.0 * M_PI), (M_PI_2 - latitude) / M_PI);
}

inline double ImageField::value_at(double u, double v) const {
    double column = u * static_cast<double>(_width) - 0.5;
    double row = v * static_cast<double>(_height) - 0.5;

    long row_low = static_cast<long>(std::floor(row));
    long column_low = static_cast<long>(std::floor(column));
    double row_fraction = row - static_cast<double>(row_low);
    double column_fraction = column - static_cast<double>(column_low);

    double brightness =
        (1.0 - row_fraction) * (
            (1.0 - column_fraction) * brightness_at(row_low, column_low) +
            column_fraction * brightness_at(row_low, column_low + 1)
        ) +
        row_fraction * (
            (1.0 - column_fraction) * brightness_at(row_low + 1, column_low) +
            column_fraction * brightness_at(row_low + 1, column_low + 1)
        );

    return _density_range.low() + (1.0 - brightness) * _density_range.measure();
}

inline double ImageField::brightness_at(long row, long column) const {
    long clamped_row = std::clamp(row, 0L, static_cast<long>(_height) - 1);
    long wrapped_column = ((column % static_cast<long>(_width)) + static_cast<long>(_width)) %
        static_cast<long>(_width);

    return _brightness[static_cast<size_t>(clamped_row) * _width + static_cast<size_t>(wrapped_column)];
}

} // namespace geometry_art::fields::scalar

#endif //GEOMETRY_ART_FIELDS_SCALAR_IMAGE_FIELD_HPP_
