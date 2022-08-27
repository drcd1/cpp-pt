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
        lps.pdf *= 1.0/float(lights.size());
        return lps;
    }

    //TODO: pdf of sample
    virtual float pdf(LightSample) const {return 1.0;}

    void add(std::shared_ptr<Light> light){
        lights.push_back(light);
    }

    void build(){

    }


};


}
#endif