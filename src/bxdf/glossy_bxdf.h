#ifndef CPPPT_GLOSSY_BXDF
#define CPPPT_GLOSSY_BXDF

#include <math/math.h>
#include <math/sampler.h>
#include <shape/intersection.h>
#include <texture/texture.h>
#include <math/math.h>
#include <math/microfacet_distribution.h>

namespace cpppt{


class GlossyBxDF: public BxDF{
    private:
        float alpha;



    public:
        GlossyBxDF(float roughness): alpha(roughness2Alpha(roughness)) {}

        Vec3 eval(const Vec3& wo, const Vec3& wi,const Intersection& it) {
            Mat3 coords(it.tangent,it.bitangent,it.normal);
            Vec3 wo_loc = normalized(coords.transpose()*(wo));
            Vec3 wi_loc = normalized(coords.transpose()*(wi));
            Vec3 wh = normalized(wo_loc+wi_loc);
            MicrofacetDistribution md(alpha,alpha);
            return md.g(wo_loc,wi_loc)*md.d(wh)/(4.0*dot(wo,it.normal));

        }
        float sample(Sampler& sampler, const Vec3& wo, const Intersection& it, Vec3* sample_direction) {

            Mat3 coords(it.tangent,it.bitangent,it.normal);
            MicrofacetDistribution md(alpha,alpha);
            Vec3 wo_loc = coords.transpose()*(wo);
            Vec3 wh_loc = normalized(md.sample(wo_loc,sampler));
            Vec3 wh = normalized(coords*(wh_loc));
            float p = md.pdf(wo_loc,wh_loc);
            *sample_direction = reflect(wo,wh);
            return pdf(wo,*sample_direction,it);
        };

        float pdf(const Vec3& wo, const Vec3& wi, const Intersection& it){
            Mat3 coords(it.tangent,it.bitangent,it.normal);
            Vec3 wo_loc = coords.transpose()*(wo);

            Vec3 wh = normalized(coords.transpose()*(wo+wi));
            return MicrofacetDistribution(alpha,alpha).pdf(wo_loc, wh);
        }
        Vec3 emit(const Vec3& wo, const Intersection& it) {
            return Vec3(0.0);
        }

        static float roughness2Alpha(float roughness) {
            return roughness*roughness;
        }
/*
        bool is_delta() override {
            return true;
        }
*/


};
}
#endif