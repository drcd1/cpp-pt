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
    bool delta;
    bool infinite;
    const Light* ref;
};

struct LightPathStart
{
    float pdf;
    float area_pdf;
    float angle_pdf;
    Vec3 position;
    Vec3 direction;
    Vec3 radiance;
};

class Light{
private:
    int id = -1;
    Light* parent = nullptr;

public:

    virtual LightSample connect_eye_path(Sampler&, const Intersection&) const= 0;
    //virtual Ray generate_light_ray(Sampler&);

    virtual LightPathStart sample(Sampler&) const = 0;


    virtual bool is_delta() const {
        return true;
    }
    virtual float pdf(int l_id, const Light* parent, const Vec3& coords, const Vec3& lit) const {
        return 1.0;
    }

    //only for infinite non delta lights
    virtual Vec3 emit(const Vec3& dir) const{
        throw std::runtime_error("Not implemented");
    }

    virtual bool is_infinite_not_delta() const{
        return false;
    }

    void set_group(int a_id, Light* a_parent){
        id = a_id;
        parent = a_parent;
    }

    std::pair<int,const Light*> get_group() const {
        return {id,parent};
    }

};


}
#endif