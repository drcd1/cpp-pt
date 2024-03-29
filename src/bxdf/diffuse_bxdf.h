#ifndef CPPPT_DIFFUSE_BXDF
#define CPPPT_DIFFUSE_BXDF

#include <math/math.h>
#include <math/sampler.h>
#include <shape/intersection.h>
#include <texture/texture.h>
#include <math/math.h>

namespace cpppt{


class DiffuseBxDF: public BxDF{

    private:
        Vec3 col;
    public:

        DiffuseBxDF(Vec3 col):col(col){
            #ifdef _DEBUG
            name = "DiffuseBxDF";
            #endif
        }

        Vec3 eval(const Vec3& wo, const Vec3& wi, const Intersection& it) {
            return col*fabs((dot(it.normal, wi)))/M_PI;
        }
        DirectionalSample sample(Sampler& sampler, const Vec3& wo, const Intersection& it) {
            float r1 = sampler.sample();
            float r2 = sampler.sample();
            Vec3 sample = sample_hemisphere_cos(r1, r2);
            float p = sample.z/M_PI;

            Vec3 n = correct_normal(it.normal,wo);

            Vec3 x,y,z;
            orthogonal(n,&x,&y,&z);
            Mat3 coords(y,z,x);
            Vec3 sample_direction = coords*sample;
            return DirectionalSample(sample_direction,p,false);

        };
        Vec3 emit(const Vec3& wo, const Intersection& it) {
            return Vec3(0.0);
        }

        float pdf(const Vec3& wo, const Vec3& wi,const Intersection& it) {
            return fabs(dot(it.normal,wi))/M_PI;
        }

};
}
#endif