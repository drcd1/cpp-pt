#ifndef CPPPT_MIRROR_BXDF
#define CPPPT_MIRROR_BXDF

#include <math/math.h>
#include <math/sampler.h>
#include <shape/intersection.h>
#include <texture/texture.h>
#include <math/math.h>

namespace cpppt{


class MirrorBxDF: public BxDF{

    public:
        Vec3 eval(const Vec3& wo, const Vec3& wi,const Intersection& it) {
            return Vec3(1.0,1.0,1.0);
        }
        float sample(Sampler& sampler, const Vec3& wo, const Intersection& it, Vec3* sample_direction) {

            *sample_direction = reflect(wo,it.normal);
            return 1.0;
        };

        float pdf(const Vec3& wo, const Vec3& wi,const Intersection& it) {
            return 1.0;
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