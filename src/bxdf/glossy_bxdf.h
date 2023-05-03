#ifndef CPPPT_GLOSSY_BXDF
#define CPPPT_GLOSSY_BXDF

#include <math/math.h>
#include <math/sampler.h>
#include <shape/intersection.h>
#include <texture/texture.h>
#include <math/math.h>
#include <math/microfacet_distribution.h>
#include <algorithm>

namespace cpppt{


class GlossyBxDF: public BxDF{
    private:
        float alpha;
        bool delta;



    public:
        GlossyBxDF(float roughness): alpha(roughness2Alpha(roughness)), delta(roughness<0.01) {
            #ifdef _DEBUG
            name = "GlossyBxDF";
            #endif
        }

        Vec3 eval(const Vec3& wo, const Vec3& wi,const Intersection& it) const override {
            Mat3 coords(it.tangent,it.bitangent,it.normal);
            Vec3 wo_loc = normalized(coords.transpose()*(wo));
            Vec3 wi_loc = normalized(coords.transpose()*(wi));
            Vec3 wh = normalized(wo_loc+wi_loc);
            MicrofacetDistribution md(alpha,alpha);
            float e = md.g(wo_loc,wi_loc)*md.d(wh)/(4.0*dot(wo,it.normal));
            if(std::isnan(e))
                std::cout<<"nannn21"<<std::endl;
            
            if(e<0.0f)
                e = 0.0f;
            return Vec3(e);

        }
        DirectionalSample sample(Sampler& sampler, const Vec3& wo, const Intersection& it) const override {
            if(delta){
                return DirectionalSample(reflect(wo,it.normal),1.0,true);
            }
            Mat3 coords(it.tangent,it.bitangent,it.normal);
            MicrofacetDistribution md(alpha,alpha);
            Vec3 wo_loc = normalized(coords.transpose()*(wo));
            Vec3 wh_loc = normalized(md.sample(wo_loc,sampler));
            Vec3 wh = normalized(coords*(wh_loc));
            float p = md.pdf(wo_loc,wh_loc);
            Vec3 sample_direction = reflect(wo,wh);
            return DirectionalSample(sample_direction, pdf(wo,sample_direction,it), false);
        };

        float pdf(const Vec3& wo, const Vec3& wi, const Intersection& it) const override {
            Mat3 coords(it.tangent,it.bitangent,it.normal);
            //todo: check if matrix is orthonormal
            Vec3 wo_loc = normalized(coords.transpose()*(wo));
            Vec3 wh = normalized(coords.transpose()*(wo+wi));
            /*if(dot(wo,it.normal)*dot(wi,it.normal)<=0.0){
                //todo: check this?
                return 0.0;
            }*/

            float p = MicrofacetDistribution(alpha,alpha).pdf(wo_loc, wh)/(4.0*fabs(dot(wo_loc,wh)));

            if(std::isnan(p)){
                std::cout<<"nan again";
            }

            return p;
        }
        Vec3 emit(const Vec3& wo, const Intersection& it) const override {
            return Vec3(0.0);
        }

        static float roughness2Alpha(float roughness) {
            return std::max(roughness*roughness,0.01f);
        }
/*
        bool is_delta() override {
            return true;
        }
*/


};
}
#endif