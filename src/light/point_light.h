#ifndef CPPPT_POINT_LIGHT
#define CPPPT_POINT_LIGHT

#include <light/light.h>

namespace cpppt{
class PointLight: public Light {
private:
    Vec3 position;
    Vec3 color;
    float intensity;

public:
    PointLight(Vec3 position, Vec3 color, float intensity): position(position),color(color),intensity(intensity){}
    //for NEE
    LightSample connect_eye_path(Sampler& s, const Intersection& other) const {
        LightSample ls;
        ls.position = position;
        ls.pdf = 1.0;
        ls.intensity = color*intensity/(4.0*M_PI);
        ls.normal = Vec3(1.0,0.0,0.0);
        ls.ref = this;
        return ls;
        /*TODO: divide by 4pi??*/
    }

    LightPathStart sample(Sampler& s) const {
        LightPathStart lps;
        lps.pdf = 1.0/(4.0*M_PI);
        lps.radiance = color*intensity; //TODO: divide 4pi?
        lps.position = position;
        lps.direction = sample_sphere(s.sample(),s.sample());
        return lps;
    }
    /*
    virtual Ray generate_light_ray(Sampler& s) {
        float r1 = s.sample();
        float r2 = s.sample();
        throw "not implemented";
    }*/
};
}
#endif