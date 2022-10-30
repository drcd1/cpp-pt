#ifndef CPPPT_DISNEY_BXDF
#define CPPPT_DISNEY_BXDF

#include <math/math.h>
#include <math/sampler.h>
#include <shape/intersection.h>
#include <texture/texture.h>
#include <math/math.h>
#include <bxdf/glossy_bxdf.h>

namespace cpppt{

float fresnel_schlick(float r0, float cos_h){
    if(r0 == 0.0){
        return 0.0;
    }
    float tmp = 1.0-cos_h;
    float tmp2 = tmp*tmp;
    return r0 + (1.0-r0)*tmp2*tmp2*tmp;
}

Vec3 compute_halfway_vector(Vec3 wo, Vec3 wi){
    return normalized(wo+wi);
}

float ior2spec(float ior){
    float a = ior-1;
    float b = ior+1;
    a = a/b;
    a = a*a;
    return a;
}


/*note: this spec value has been divided by 0.08 in standard material*/
class DisneyBxDF: public BxDF{
    private:
        Vec3 albedo;
        float metalness;
        float roughness;
        float specular;
        float transparent;
        float ior;
        Vec3 transparent_tint;
        bool delta;

        //
        //Cook-torrance: G()*D()/(cos(wi)*cos(wo))
    public:
        DisneyBxDF(Vec3 albedo, float metalness, float roughness, float specular, float transparent,float ior):
        albedo(albedo), metalness(metalness),
        roughness(roughness), specular(specular*0.08), transparent(transparent),ior(ior), delta(roughness<0.01?true:false) {}


        Vec3 eval(const Vec3& wo, const Vec3& wi,const Intersection& it) {

            auto it2 = it;
            if(dot(it.g_normal, wo)<0.0){
                it2.g_normal = it.g_normal*(-1.0);
                it2.normal = it.normal*(-1.0);
            }
            Vec3 ret(0.0);
            if( same_hemisphere(it2,wi,wo)){
                Vec3 hv = compute_halfway_vector(wi,wo);
                float fr0 = fresnel_schlick(specular, dot(hv,wi));
                Vec3 fr0_m;

                float m = metalness;
                float s = (1.0-metalness)*(fr0);
                float t = transparent * (1.0-fr0) * (1.0-metalness);

                if(m>0.0){
                    fr0_m.x = fresnel_schlick(albedo.x, dot(hv,wi));
                    fr0_m.y = fresnel_schlick(albedo.y, dot(hv,wi));
                    fr0_m.z = fresnel_schlick(albedo.z, dot(hv,wi));
                }


                /*TODO: weight diffuse better*/
                float d = (1.0-(fr0)) * (1.0-metalness) * (1.0 -transparent);


                //Todo: check s + t + m + d = 1.0
                if(s>EPS || m>EPS){
                    Vec3 tmp(1.0);
                    if(!delta){
                        tmp = GlossyBxDF(roughness).eval(wo,wi,it2);
                    }
                    ret = ret + tmp*s + tmp*fr0_m*m;
                }

                if(d>EPS){
                    ret = ret + DiffuseBxDF(albedo).eval(wo,wi,it2)*d;
                }

            }
           // if(t>EPS){
                /*
                if(roughness>0.0)
                    ret += transparent_tint*t*TransparentGlossyBxDF(roughness).eval(wo,wi,it);
                //else don't eval specular: returns zero;
                */
           // }


            return ret;

        }
        BxDFSample sample(Sampler& sampler, const Vec3& wo, const Intersection& it) {
            Vec3 diff_sample;
            Vec3 spec_sample;
            Vec3 transp_sample;
            auto it2 = it;
            if(dot(it.g_normal, wo)<0.0){
                it2.g_normal = it.g_normal*(-1.0);
                it2.normal = it.normal*(-1.0);
            }

            DiffuseBxDF diff(albedo);
            GlossyBxDF glossy(roughness);
            Vec3 sd_g;
            BxDFSample gl_samp = glossy.sample(sampler,wo,it2);
            sd_g = gl_samp.wi;
            float pg = gl_samp.pdf;
            Vec3 hv = compute_halfway_vector(sd_g,wo);
            float fr0 = fresnel_schlick(specular, dot(hv,sd_g));


            float s = fr0;
            float t = transparent * (1.0-fr0);
            float m = metalness * (1.0-transparent) * (1.0-fr0);

            /*TODO: weight diffuse better*/
            float d = (1.0-(fr0)) * (1.0-metalness) * (1.0 -transparent);
            d = 1.0;
            //just choose a random based on weights
            float r = sampler.sample();
            if(r<(s+m)){ //reflect
                return BxDFSample(sd_g,(s+m)*pg,gl_samp.delta);
            } else if(r<s+m+t){ //refract
            /*
                Intersection new_inter = intersection;
                new_inter.normal = hv;
                //todo:check if prob is the same for refraction
                return RefractionBxDF(ior).sample(sampler,wo,Intersection,sample_direction)*pg*t
            */
                //RefractionBxDF rbxdf(ior);
                //BxDFSample refraction_sample = rbxdf.sample(sampler, wo,it);
                return BxDFSample(sd_g, 1.0,true);
                //return refraction_sample;

            } else {

                BxDFSample diff_samp =  diff.sample(sampler,wo,it2);
                diff_samp.pdf *= d;
                return diff_samp;
            }

        };

        float pdf(const Vec3& wo, const Vec3& wi, const Intersection& it){
            //if wi is on the other side, compute the hv and wi as reflect

            auto it2 = it;
            if(dot(it.g_normal, wo)<0.0){
                it2.g_normal = it.g_normal*(-1.0);
                it2.normal = it.normal*(-1.0);
            }

            if(!non_zero(it2, wi,wo)){
                return 0.0;
            }
            Vec3 hv = normalized(wi+wo);
            float fr0 = fresnel_schlick(specular, dot(hv,wi));

            float s = fr0;
            float t = transparent * (1.0-fr0);
            float m = metalness * (1.0-transparent) * (1.0-fr0);
            //todo: weight diffuse better
            float d = (1.0-(fr0)) * (1.0-metalness) * (1.0 -transparent);
            float ret = DiffuseBxDF(albedo).pdf(wo,wi,it2)*d +  GlossyBxDF(roughness).pdf(wo,wi,it2)*(s+m);

            return ret;

        }

        virtual Vec3 emit(const Vec3& wo, const Intersection& it) {
            return Vec3(0.0);
        }

};
}
#endif