#ifndef CPPPT_MIRROR_BXDF
#define CPPPT_MIRROR_BXDF

#include <math/math.h>
#include <math/sampler.h>
#include <shape/intersection.h>
#include <texture/texture.h>
#include <math/math.h>
#include <bxdf/bxdf.h>

namespace cpppt{


class MirrorBxDF: public BxDF{

    public:
        Vec3 eval(const Vec3& wo, const Vec3& wi,const Intersection& it) const override {
            return Vec3(1.0,1.0,1.0);
        }
        DirectionalSample sample(Sampler& sampler, const Vec3& wo, const Intersection& it) const override {
            DirectionalSample tmp;
            tmp.wi =   reflect(wo,it.normal);
            tmp.delta = true;
            tmp.pdf = 1.0;
            return tmp;
        };

        float pdf(const Vec3& wo, const Vec3& wi,const Intersection& it) const override {
            return 1.0;
        }

        Vec3 emit(const Vec3& wo, const Intersection& it) const override {
            return Vec3(0.0);
        }


};
}
#endif