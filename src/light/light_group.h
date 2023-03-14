#ifndef CPPPT_LIGHT_GROUP
#define CPPPT_LIGHT_GROUP

#include <math/math.h>
#include <math/sampler.h>
#include <shape/intersection.h>
#include <texture/texture.h>
#include <math/math.h>
#include <primitive/ray.h>
#include <light/light.h>
#include <light/environment_light.h>

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


    LightPathStart sample(Sampler& s) const {
        int sampled = int(s.sample()*lights.size())%lights.size();

        LightPathStart lps = lights.at(sampled)->sample(s);
        lps.radiance = lps.radiance;
        lps.pdf *= 1.0/float(lights.size());
        lps.area_pdf *= 1.0/float(lights.size());
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

    int size(){
        return lights.size();
    }

    void build(){

    }


};

class SceneLights: public Light{
private:
    std::shared_ptr<LightGroup> lg;
    std::shared_ptr<EnvironmentLight> env;
public:
    SceneLights(std::shared_ptr<LightGroup> lg, std::shared_ptr<EnvironmentLight> env):lg(lg),env(env){
    }
    virtual LightSample connect_eye_path(Sampler& s, const Intersection& hit) const{
        float r = s.sample();
        LightSample ls;
        if(r>0.5){
            ls = lg->connect_eye_path(s,hit);
        } else {
            ls = env->connect_eye_path(s,hit);
        }
        ls.pdf *= 0.5;
        return ls;
    };


    LightPathStart sample(Sampler& s) const {
        float r = s.sample();
        LightPathStart lps;
        if(r>0.5){
            lps = lg->sample(s);
        } else {
            lps = env->sample(s);
        }
        lps.pdf*=0.5;
        return lps;
    }

    float pdf(int l_id, const Light* parent, const Vec3& coords, const Vec3& lit) const {
        if(l_id<0){
            return 0.5*env->pdf(l_id,parent,coords,lit);
        } else {
            return 0.5*lg->pdf(l_id, parent, coords, lit);
        }
    }

    Vec3 emit(const Vec3& dir)const override {
        return env->emit(dir);
    }

    bool is_infinite_not_delta()const override {
        return true;
    }
};


}
#endif