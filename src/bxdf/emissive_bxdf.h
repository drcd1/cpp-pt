#ifndef CPPPT_EMISSIVE_BXDF
#define CPPPT_EMISSIVE_BXDF

#include <math/math.h>
#include <math/sampler.h>
#include <shape/intersection.h>
#include <texture/texture.h>

namespace cpppt{


class EmissiveBxDF: public BxDF{
    private:
        std::shared_ptr<Texture<Vec3>> tex;
        float intensity;
    public:
        EmissiveBxDF(std::shared_ptr<Texture<Vec3>> tex, float intensity): tex(tex),intensity(intensity){}

        Vec3 eval(const Vec3& wo, const Vec3& wi,const Intersection& it) {
            return Vec3(0.0,0.0,0.0);
        }
        float sample(Sampler& sampler, const Vec3& wo, const Intersection& it, Vec3* sample_direction) {
            return -1.0;
        };
        Vec3 emit(const Vec3& wo, const Intersection& it) {
            return tex->sample(it.texture_coords)*intensity;
        }

};
}
#endif