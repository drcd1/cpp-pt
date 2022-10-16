#ifndef CPPPT_LIGHT_GROUP
#define CPPPT_LIGHT_GROUP

#include <math/math.h>
#include <math/sampler.h>
#include <shape/intersection.h>
#include <texture/texture.h>
#include <math/math.h>
#include <primitive/ray.h>
#include <light/light.h>

namespace cpppt{



class LightGroup: public Light{
private:
    std::vector<std::shared_ptr<Light>> lights;
public:
    /*
    virtual LightSample sample(Sampler& s, const Intersection& hit) const{
        //uniformly sample

    };*/
    virtual LightSample connect_eye_path(Sampler& s, const Intersection& hit) const{
        int sampled = int(s.sample()*lights.size())%lights.size();
        LightSample ls = lights.at(sampled)->connect_eye_path(s,hit);
        ls.pdf *= 1.0/float(lights.size());
        return ls;
    };
  //  virtual Ray generate_light_ray(Sampler&) = 0;

    LightPathStart sample(Sampler& s) const {
        int sampled = int(s.sample()*lights.size())%lights.size();

        LightPathStart lps = lights.at(sampled)->sample(s);
        lps.radiance = lps.radiance*lights.size();
        lps.pdf *= 1.0/float(lights.size());
        return lps;
    }

    //TODO: pdf of sample
    float pdf(int l_id, const Light* parent, const Vec3& coords, const Vec3& lit) const {
        if(parent!=this){
            throw std::runtime_error("PDF of light group from light in wrong group. TODO: fix");
        }
        auto light = lights.at(l_id);
        float pdf = 1.0/float(lights.size());
        if(light->is_delta()){
            return pdf;
        }

        pdf *= light->pdf(l_id, parent, coords, lit);
        return pdf;
    }

    void add(std::shared_ptr<Light> light){
        light->set_group(lights.size(), this);
        lights.push_back(light);
    }

    void build(){

    }


};


}
#endif