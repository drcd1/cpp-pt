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


        virtual bool non_zero(const Intersection& i, const Vec3& wi, const Vec3& wo){
            return same_hemisphere(i,wi,wo);
        }

        Vec3 eval(const Vec3& wo, const Vec3& wi,const Intersection& it) {

            Vec3 ret(0.0);
            if( same_hemisphere(it,wi,wo)){
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
                        tmp = GlossyBxDF(roughness).eval(wo,wi,it);
                    }
                    ret = ret + tmp*s + tmp*fr0_m*m;
                }

                if(d>EPS){
                    ret = ret + DiffuseBxDF(albedo).eval(wo,wi,it)*d;
                }
                if(std::isnan(ret.x) || std::isnan(ret.y) || std::isnan(ret.z) ){
                    std::cout<<"NAN3!"<<std::endl;
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
        virtual float sample(Sampler& sampler, const Vec3& wo, const Intersection& it, Vec3* sample_direction) {
            Vec3 diff_sample;
            Vec3 spec_sample;
            Vec3 transp_sample;

            DiffuseBxDF diff(albedo);
            GlossyBxDF glossy(roughness);
            Vec3 sd_g;
            float pg = glossy.sample(sampler,wo,it,&sd_g);
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
                *sample_direction = sd_g;
                if(std::isnan((s+m)*pg)){
                    std::cout<<"NAN53"<<std::endl;
                }
                return (s+m)*pg;
            } else if(r<s+m+t){ //refract
            /*
                Intersection new_inter = intersection;
                new_inter.normal = hv;
                //todo:check if prob is the same for refraction
                return RefractionBxDF(ior).sample(sampler,wo,Intersection,sample_direction)*pg*t
            */

            return 0.0;

            } else {
                float p =  diff.sample(sampler,wo,it,sample_direction)*d;
                if(std::isnan(p)){
                    std::cout<<"NAN523"<<std::endl;
                }

                return p;
            }

        };

        float pdf(const Vec3& wo, const Vec3& wi, const Intersection& it){
            //if wi is on the other side, compute the hv and wi as reflect
            if(non_zero(it, wi,wo)){
                return 0.0;
            }
            Vec3 hv = wi+wo;
            float fr0 = fresnel_schlick(specular, dot(hv,wi));

            float s = fr0;
            float t = transparent * (1.0-fr0);
            float m = metalness * (1.0-transparent) * (1.0-fr0);
            //todo: weight diffuse better
            float d = (1.0-(fr0)) * (1.0-metalness) * (1.0 -transparent);
            float ret = DiffuseBxDF(albedo).pdf(wo,wi,it)*d +  GlossyBxDF(roughness).pdf(wo,wi,it)*(s+m);
            if(std::isnan(ret)){
                std::cout<<"NAN4!<"<<std::endl;
            }

            return ret;

        }

        virtual Vec3 emit(const Vec3& wo, const Intersection& it) {
            return Vec3(0.0);
        }

};
}
#endif