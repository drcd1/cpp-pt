#ifndef CPPPT_PATHTRACER_MIS
#define CPPPT_PATHTRACER_MIS

#include <renderer/renderer.h>
#include <math/math.h>
#include <camera/camera.h>
#include <image/rgb_image.h>
#include <primitive/ray.h>
#include <shape/intersection.h>
#include <primitive/primitive_leaf.h>
#include <scene.h>
#include <math/sampler.h>
#include <omp.h>
#include <iostream>
#include <light/shape_light.h>
namespace cpppt{
class PathtracerMIS : public Renderer{
    int samples;

    Vec3 render_sky(const Ray& r) const {
        return Vec3(.0);
    }

    static float russian_roulette(const Vec3& col){
        //Todo: better rr
        return (fabs(col.x*0.2 + col.y*0.5 +col.z*0.3))*0.5 + 0.4;
    }

    Vec3 integrate(const Scene& scene, const Vec2& coords, Sampler& sampler) const {
        Ray ray = scene.camera->get_ray(coords);
        Intersection intersection;
        Vec3 col(0.0);
        Vec3 mul(1.0);
        bool sampled_delta = true;
        float p = 1.0;
        for(int i = 0; i<32; i++){
            bool intersected = scene.primitive->intersect(ray,&intersection);
            if(!intersected){
                col = col + mul*render_sky(ray);
                break;
            } else {
                auto bsdf = intersection.get_bxdf();
                if(bsdf->is_emitter()){
                    if(sampled_delta){

                        //scene.light = dynamic_cast
                        //float p_bsdf_nee = scene.light->pdf(static_cast<std::shared_ptr<EmissiveBxDF>>(bsdf))
                        //float p_bsdf_nee  = scene.light->pdf()*(r*r)/bsdf_cos;

                        Vec3 rad = bsdf->emit(ray.d*(-1.0),intersection);

                        col = col + mul*rad;
                    } else {
                        float p_bsdf_bsdf = p;

                        auto lid = intersection.primitive->get_light()->get_group();

                        float r = ray.max_t;
                        float cos_theta = fabs(dot(intersection.normal,ray.d)) + EPS;
                        float p_bsdf_nee = scene.light->pdf(lid.first,lid.second, intersection.texture_coords,ray.o)*r*r/cos_theta;

                        float w_mis = p_bsdf_bsdf/(p_bsdf_bsdf+p_bsdf_nee);

                        Vec3 rad = bsdf->emit(ray.d*(-1.0),intersection);
                        col = col + mul*rad*w_mis;
                    }
                    break; //emitters do not reflect!
                }
                sampled_delta = false;



                Vec3 sample_direction;
                //float w_bsdf;

                p = bsdf->sample(sampler, ray.d*(-1.0), intersection, &sample_direction);

                //TODO: if sample is not valid, full absorption

                /* NEE */
                if(!bsdf->is_delta()){
                    LightSample light_sample = scene.light->connect_eye_path(sampler, intersection);
                    Vec3 s_dir = (light_sample.position-intersection.hitpoint);
                    float len = length(s_dir);
                    s_dir = s_dir/len;

                    Ray shadow_ray(intersection.hitpoint + s_dir*EPS,
                            s_dir,
                            len-2.0*EPS);


                    if(bsdf->non_zero(intersection,ray.d*(-1.0),shadow_ray.d)){
                    if(!scene.primitive->intersect_any(shadow_ray)){


                        //float cosine_term = 1.0;
                        //if(!light_sample.ref->is_delta())
                        //    float cosine_term = fabs(dot(light_sample.normal,shadow_ray.d));
                       /*TODO:why doesn't the cosine term here work?*/
                        float r = len;
                        float p_nee_bsdf = bsdf->pdf(ray.d*(-1.0),s_dir,intersection);
                        float cosine_term = fabs(dot(light_sample.normal,s_dir))+EPS;
                        float p_nee_nee = light_sample.pdf*r*r/cosine_term;
                        float w_mis = p_nee_nee/(p_nee_bsdf + p_nee_nee);

                        Vec3 eval = bsdf->eval(
                            ray.d*(-1.0),
                            shadow_ray.d,
                            intersection
                        );
                        col = col + mul*eval*light_sample.intensity/(light_sample.pdf*r*r)*w_mis;
                    }
                    }
                } else {
                    sampled_delta = true;
                }
                /* MIS */

                //float p_nee_bsdf = bsdf->pdf(intersection, light_sample);
                //float p_nee_nee = light_sample.p*r*r/some_cosine;

                //float w_bsdf = p_nee_nee/(p_nee_bsdf + p_nee_nee);

                if(!bsdf->non_zero(intersection,ray.d*(-1.0),sample_direction)){
                    break;
                }

                Vec3 eval = bsdf->eval(ray.d*(-1.0), sample_direction, intersection);

                if(i>2){
                    float rr = russian_roulette(eval);
                    if(sampler.sample() > rr) {
                        break;
                    }
                    p = p*rr;
                }
                mul = mul*eval/p;
                ray = Ray(intersection.hitpoint+sample_direction*EPS, sample_direction);
            }
        }
        return col;
    }

public:
    PathtracerMIS(const RenderSettings& rs): samples(rs.spp) {}

    void render(Scene& sc, std::string filename) const {
        RgbImage* image= &(sc.camera->get_image());
        Vec2i res = image->res;

        #pragma omp parallel for
        for(int i = 0; i<res.x; i++){


            RandomSampler s(i);
            if(i%50==0)
                std::cout<<"rendering line "<<i<<std::endl;
            for(int j = 0; j<res.y; j++){
                Vec3 acc(0.0);
                for(int k = 0; k<samples; k++){
                    for(int l = 0; l<samples; l++){
                        float r1 = s.sample();
                        float r2 = s.sample();

                        Vector2<float> coords( ((float(i)+float(k+r1)/samples)/float(res.x))*2.0-1.0, -(((float(j)+float(l+r2)/samples)/float(res.y))*2.0-1.0));
                        acc = acc +integrate(sc, coords, s);

                    }
                }

                acc = acc/(float(samples*samples));
                image->put_pixel(i,j,acc);

            }
        }


        image->save(filename);
    }
};
}
#endif