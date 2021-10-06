#ifndef CPPPT_CAMERA
#define CPPPT_CAMERA

#include <math/math.h>
#include <primitive/ray.h>
#include <image/rgb_image.h>
namespace cpppt{

class Camera {
public:
    virtual Ray get_ray(Vec2 coords) const = 0;
    virtual RgbImage& get_image() = 0;
};

}


#endif