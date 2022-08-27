#ifndef CPPPT_LIGHT
#define CPPPT_LIGHT

#include <math/math.h>
#include <math/sampler.h>
#include <shape/intersection.h>
#include <texture/texture.h>
#include <math/math.h>
#include <primitive/ray.h>

namespace cpppt{


class Light;

struct LightSample
{
    float pdf;
    Vec3 position;
    Vec3 intensity;
    Vec3 normal;
    const Light* ref;
};

struct LightPathStart
{
    float pdf;
    Vec3 position;
    Vec3 direction;
    Vec3 radiance;
};

class Light{

public:

    virtual LightSample connect_eye_path(Sampler&, const Intersection&) const= 0;
    //virtual Ray generate_light_ray(Sampler&);


    //pdf of
    virtual float pdf(LightSample) const {return 1.0;}

    virtual LightPathStart sample(Sampler&) const = 0;


    virtual bool is_delta() const {
        return true;
    }

};


}
#endif