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

    #ifndef NO_GUI
    std::vector<int> counter;
    int total_pixels;
    int current_pixels;
    #endif

     Vec3 render_sky(const Ray& r, const Scene& s) const {
        if(s.light->is_infinite_not_delta()){
            return s.light->emit(r.d);
        } else {
            return Vec3(.0);
        }
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
                Vec3 ret = render_sky(ray,scene);
                float p_bsdf_bsdf = p;
                float p_bsdf_nee = scene.light->pdf(-1,scene.light.get(), ray.d,ray.o);
                float mis = p_bsdf_bsdf/(p_bsdf_nee+p_bsdf_bsdf);
                //TODO: NO MIS WEIGHT IF i = 0
                if(i==0){
                    mis = 1.0;
                }
                col = col + mul*ret*mis;
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



                DirectionalSample sample = bsdf->sample(sampler, ray.d*(-1.0), intersection);

                Vec3 sample_direction = sample.wi;
                p = sample.pdf;


                /* NEE */
                if(!sample.delta){
                    LightSample light_sample = scene.light->connect_eye_path(sampler, intersection);
                    Vec3 s_dir = (light_sample.position-intersection.hitpoint);
                    float len = length(s_dir);
                    s_dir = s_dir/len;
                    if(light_sample.infinite){
                        len = 10e6;
                    }

                    Ray shadow_ray(intersection.hitpoint + s_dir*EPS,
                            s_dir,
                            len-2.0*EPS);


                    if(bsdf->non_zero(intersection,ray.d*(-1.0),shadow_ray.d)){
                    if(!scene.primitive->intersect_any(shadow_ray)){
                        float r = len;
                        if(light_sample.infinite){
                            r = 1.0;
                        }
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

#ifndef NO_GUI
    void updateProgress(RenderProgress& progress) override{
        total_pixels = std::max(total_pixels,1);
        progress.percentage = float(current_pixels)/total_pixels;
        //progress.percentage= 0.3;
        progress.current_task = 1;
        progress.total_tasks = 1;
        sprintf(progress.task_description,"Rendering");
    }
#endif

    PathtracerMIS(const RenderSettings& rs): samples(rs.spp) {}

    void render(Scene& sc, std::string filename) {
        RgbImage* image= &(sc.camera->get_image());
        Vec2i res = image->res;

        #ifndef NO_GUI
        int counter_size = omp_get_max_threads();
        counter.resize(counter_size);
        total_pixels = res.x*res.y;
        #endif

        #pragma omp parallel for
        for(int i = 0; i<res.x; i++){
            counter.at(omp_get_thread_num()) += res.y;


            RandomSampler s(i);

            #ifndef NO_GUI

            if(omp_get_thread_num()==0){
                int tmp = 0;
                for(int k = 0; k<counter_size; k++){
                    tmp+=counter.at(k);
                }
                current_pixels = tmp;
            }
            #endif

            //if(i%50==0)
            //    std::cout<<"rendering line "<<i<<std::endl;
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
    static const char* name(){
        return "PathtracerMIS";
    }
};
}
#endif