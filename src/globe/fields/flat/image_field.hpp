#ifndef GLOBEART_SRC_GLOBE_FIELDS_FLAT_IMAGE_FIELD_HPP_
#define GLOBEART_SRC_GLOBE_FIELDS_FLAT_IMAGE_FIELD_HPP_

#include "../scalar/image_field.hpp"
#include "../../math/interval.hpp"
#include "../../types.hpp"
#include <string>
#include <utility>

namespace globe::fields::flat {

using globe::math::Interval;

// An image draped over the rectangle: x runs across the picture, y runs up
// it, darker is denser. Periodicity across the seams is the image's own
// affair -- a tileable image wraps invisibly, any other shows its cut.
class ImageField {
 public:
    ImageField(scalar::ImageField image, double width, double height);

    [[nodiscard]] static ImageField load(
        const std::string& path,
        Interval density_range,
        double width,
        double height
    );

    [[nodiscard]] double value(const Vector2& point) const;

 private:
    scalar::ImageField _image;
    double _width;
    double _height;
};

inline ImageField::ImageField(scalar::ImageField image, double width, double height) :
    _image(std::move(image)),
    _width(width),
    _height(height) {
}

inline ImageField ImageField::load(
    const std::string& path,
    Interval density_range,
    double width,
    double height
) {
    return ImageField(scalar::ImageField::load(path, density_range), width, height);
}

inline double ImageField::value(const Vector2& point) const {
    return _image.value_at(point.x() / _width, 1.0 - point.y() / _height);
}

} // namespace globe::fields::flat

#endif //GLOBEART_SRC_GLOBE_FIELDS_FLAT_IMAGE_FIELD_HPP_
