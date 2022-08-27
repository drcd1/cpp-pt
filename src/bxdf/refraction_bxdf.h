#ifndef CPPPT_REFRACTION_BXDF
#define CPPPT_REFRACTION_BXDF

#include <math/math.h>
#include <math/sampler.h>
#include <shape/intersection.h>
#include <texture/texture.h>
#include <math/math.h>

namespace cpppt{


class RefractionBxDF: public BxDF{

    private:
        float ior;
    public:
        RefractionBxDF(float ior):ior(ior){}
        Vec3 eval(const Vec3& wo, const Vec3& wi,const Intersection& it) {
            return Vec3(-1.0,-1.0,-1.0);
        }
        float sample(Sampler& sampler, const Vec3& wo, const Intersection& it, Vec3* sample_direction) {
            //Vec3 n = correct_normal(it.normal, wo);

            Vec3 n = it.normal;
            float n2 = 1.0/ior;
            if(dot(n,wo)<0.0){
                n2 = 1.0/n2;
                n = n*(-1.0);
            }
/*
            //from scratchapixel
            float c1  = dot(n,wo);
            float c2  = sqrt(1.0-n2*n2*(1.0-c1*c1));
            Vec3 r = normalized((-wo+n*c1)*n2-n*c2);
            //r = wo*(-1.0);
            */

            float cosThetaI = dot(wo,n);
            float cosThetaTSqr = 1.0 - n2*n2*(1.0-cosThetaI*cosThetaI);
            Vec3 r = normalized(wo*(-1.0)*n2 + n*(n2*cosThetaI - sqrt(cosThetaTSqr)));
            if(cosThetaTSqr <0){
                r = normalized(wo*(-1.0) + it.normal*dot(it.normal,wo)*2.0);
            }
            *sample_direction = r;
            //*sample_direction = normalized(wo*(-1.0));
            return -1.0;
        }

        float pdf(const Vec3& wo, const Vec3& wi,const Intersection& it) {
            return -1.0;
        }
        Vec3 emit(const Vec3& wo, const Intersection& it) {
            return Vec3(0.0);
        }
        bool is_delta(){
            return true;
        }


};
}
#endif