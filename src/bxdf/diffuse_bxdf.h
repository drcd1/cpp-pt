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
        std::shared_ptr<Texture<Vec3>> tex;
    public:
        DiffuseBxDF(std::shared_ptr<Texture<Vec3>> tex): tex(tex) {}

        Vec3 eval(const Vec3& wo, const Vec3& wi,const Intersection& it) {
            return tex->sample(it.texture_coords)*fabs(dot(it.normal, wi))/M_PI;
        }
        virtual float sample(Sampler& sampler, const Vec3& wo, const Intersection& it, Vec3* sample_direction) {
            float r1 = sampler.sample();
            float r2 = sampler.sample();
            Vec3 sample = sample_hemisphere_cos(r1, r2);
            float p = sample.z/M_PI;

            Vec3 n = it.normal;
            if (dot(n, wo) < 0.0) {
                n = n*(-1.0);
            }
            Vec3 x,y,z;
            orthogonal(n,&x,&y,&z);
            Mat3 coords(y,z,x);
            *sample_direction = coords*sample;
            return p;

        };
        virtual Vec3 emit(const Vec3& wo, const Intersection& it) {
            return Vec3(0.0);
        }

};
}
#endif