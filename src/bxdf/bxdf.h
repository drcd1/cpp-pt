#ifndef CPPPT_BXDF
#define CPPPT_BXDF
#ifndef _DEBUG
#define _DEBUG
#endif

#include <math/math.h>
#include <math/sampler.h>
#include <shape/intersection.h>
namespace cpppt{


class BxDF{
    public:

        #ifdef _DEBUG
        std::string name = "Abstract BxDF";
        #endif
        virtual Vec3 eval(const Vec3& wo, const Vec3& wi,const Intersection& it) const = 0;
        virtual DirectionalSample sample(Sampler& sampler, const Vec3& wo, const Intersection& it) const = 0;
        virtual float pdf(const Vec3& wo, const Vec3& wi,const Intersection& it) const = 0;
        virtual Vec3 emit(const Vec3& wo, const Intersection& it) const = 0;
        virtual bool is_emitter() const {return false;}

        virtual float emit_sample(Sampler& sample, const Intersection& it, Vec3* sample_direction) const {
            throw std::runtime_error("not implemented");
            return 1.0;
            float a = sample.sample();
            float b = sample.sample();
            Vec3 normal;
            normal = it.normal;
            if(a>0.5){
                a = a*2.0-1.0;
                normal = normal*(-1.0);
            } else {
                a = a*2.0;
            }


        }

        bool same_hemisphere(const Intersection& i, const Vec3& wo, const Vec3& wi) const {
            float sign = (dot(wo,i.g_normal));
            return sign*dot(wi,i.normal)>0.0 && sign*dot(wi,i.g_normal)>0.0;
        }
        //computues if they are in same hemisphere and whether that's good
        virtual bool non_zero(const Intersection& i, const Vec3& wi, const Vec3& wo) const {
            return same_hemisphere(i,wi,wo);
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