#ifndef CPPPT_CAMERA
#define CPPPT_CAMERA

#include <math/math.h>
#include <primitive/ray.h>
#include <image/rgb_image.h>
#include <math/sampler.h>
#include <shape/intersection.h>

namespace cpppt{

struct CameraConnection{
    int i;
    int j;
    Vec3 pos;
    float factor;
};

class Camera {
public:
    virtual Ray get_ray(Vec2 coords) const = 0;
    virtual Ray get_ray(Vec2 coords, float& angular_pdf) const = 0;

    virtual RgbImage& get_image() = 0;
    virtual CameraConnection connect_light_path(Sampler&, const Intersection&) const = 0;
    virtual float pdf(const Vec3& pt) const = 0;

};

}


#endif