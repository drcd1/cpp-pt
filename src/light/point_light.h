#ifndef CPPPT_POINT_LIGHT
#define CPPPT_POINT_LIGHT

#include <light/light.h>
#include <texture/texture.h>
namespace cpppt{
class PointLight: public Light {
private:
    Vec3 position;
    std::shared_ptr<Texture> color;
    float intensity;


/*
For point lights, intensity corresponds to total power.
*/
public:
    PointLight(Vec3 position, std::shared_ptr<Texture> color, float intensity): position(position),color(color),intensity(intensity){}
    //for NEE
    LightSample connect_eye_path(Sampler& s, const Intersection& other) const {
        LightSample ls;
        ls.position = position;
        ls.pdf = 1.0;
        ls.intensity = color->sample(Vec3(0.0))*intensity/(4.0*M_PI);
        ls.normal = Vec3(1.0,0.0,0.0);
        ls.ref = this;
        ls.delta = true;
        ls.infinite = false;
        return ls;
    }

    LightPathStart sample(Sampler& s) const {
        LightPathStart lps;
        lps.pdf = 1.0/(4.0*M_PI*M_PI);
        lps.angle_pdf = lps.pdf;
        lps.area_pdf = 1.0;
        lps.radiance = color->sample(Vec3(0.0))*intensity/(4.0*M_PI); //TODO: divide 4pi?
        lps.position = position;
        lps.direction = sample_sphere(s.sample(),s.sample());
        lps.light = this;
        return lps;
    }

    virtual Vec3 get_emission(const Vec3& dir, const Intersection* it = nullptr) const override {
        return  color->sample(Vec3(0.0))*intensity/(4.0*M_PI);
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