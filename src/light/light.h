#ifndef CPPPT_LIGHT
#define CPPPT_LIGHT

namespace cpppt{

class light{

public:
    struct LightSample
    {
        float pdf;
        Vec3 point;
        Vec3 intensity;
        Light* ref;
    };
    
    LightSample sample(Sampler&);
    float pdf(LightSample);

    bool specular() const {
        return true;
    }

};


}
#endif