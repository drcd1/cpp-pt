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
        virtual float pdf(const Vec3& wo, const Vec3& wi,const Intersection& it) = 0;
        virtual Vec3 emit(const Vec3& wo, const Intersection& it) = 0;
        virtual bool is_emitter(){return false;}
        virtual bool is_delta(){ return false;}

        virtual float emit_sample(Sampler& sample, const Intersection& it, Vec3* sample_direction){
            return 0.0;
        }


        Vec3 correct_normal(const Vec3& n, const Vec3& wo ) const {
            if (dot(n, wo) < 0.0) {
                return n*(-1.0);
            }
            return n;
        }


};
}
#endif