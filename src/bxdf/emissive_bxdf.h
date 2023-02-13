#ifndef CPPPT_EMISSIVE_BXDF
#define CPPPT_EMISSIVE_BXDF

#include <math/math.h>
#include <math/sampler.h>
#include <shape/intersection.h>
#include <texture/texture.h>

namespace cpppt{


class EmissiveBxDF: public BxDF{
    private:
        Vec3 col;
        bool double_sided;

    public:
        EmissiveBxDF(Vec3 col, bool double_sided=true): col(col), double_sided(double_sided){}

        Vec3 eval(const Vec3& wo, const Vec3& wi,const Intersection& it) {
            return Vec3(0.0,0.0,0.0);
        }
        DirectionalSample sample(Sampler& sampler, const Vec3& wo, const Intersection& it) {
            throw std::runtime_error("We should never do this!");
            return DirectionalSample(Vec3(0.0),-1.0,true);
        };


        float pdf(const Vec3& wo, const Vec3& wi,const Intersection& it) {
            throw std::runtime_error("We should never do this!");
            return -1.0;
        }

        Vec3 emit(const Vec3& wo, const Intersection& it) {
            return col;
        }
        virtual float emit_sample(Sampler& sampler, const Intersection& it, Vec3* sample_direction){
            float r1 = sampler.sample();
            float r2 = sampler.sample();
            Vec3 sample = sample_hemisphere_cos(r1, r2);
            float p = 0.5*sample.z/M_PI;

            //TODO: double normals?

            float r = sampler.sample();
            Vec3 n = it.normal;
            Vec3 t = it.tangent;
            if(r<0.5){
                n = n*(-1.0);
                t = t*(-1.0);
            }

            Vec3 x,y,z;
            //orthogonal(it.normal,&x,&y,&z);
            Mat3 coords(t,it.bitangent,n);
            *sample_direction = coords*sample;
            return p;
        }

        bool is_emitter(){
            return true;
        }

};
}
#endif