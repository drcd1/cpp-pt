#ifndef CPPPT_BXDF
#define CPPPT_BXDF

#include <math/math.h>
#include <math/sampler.h>
namespace cpppt{

class Intersection;

class BxDF{
    public:
        virtual Vec3 eval(const Vec3& wo, const Vec3& wi,const Intersection& it) = 0;
        virtual float sample(Sampler& sampler, const Vec3& wo, const Intersection& it, Vec3* sample_direction) = 0;
        virtual Vec3 emit(const Vec3& wo, const Intersection& it) = 0;

};
}
#endif