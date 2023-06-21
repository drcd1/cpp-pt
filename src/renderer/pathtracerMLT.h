#ifndef CPPPT_PATHTRACER_MLT
#define CPPPT_PATHTRACER_MLT

#include <renderer/renderer.h>
#include <math/math.h>
#include <camera/camera.h>
#include <image/rgb_image.h>
#include <primitive/ray.h>
#include <shape/intersection.h>
#include <scene.h>
#include <math/sampler.h>
#include <omp.h>
#include <iostream>
namespace cpppt{

class PathtracerMLT : public Renderer{
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
        Intersection prev_intersection;
        Vec3 col(0.0);
        Vec3 mul(1.0);
        bool sampled_delta = true;
        for(int i = 0; i<32; i++){
            bool intersected = scene.primitive->intersect(ray,&intersection);
            if(!intersected){
                col = col + mul*render_sky(ray);
                break;
            } else {
                auto bsdf = intersection.get_bxdf();
                if(bsdf->is_emitter()){
                    if(sampled_delta){
                        col = col + mul*bsdf->emit(ray.d*(-1.0),intersection);
                    }
                    break;
                }
                sampled_delta = false;

                DirectionalSample sample = bsdf->sample(sampler, ray.d*(-1.0), intersection);

                Vec3 sample_direction = sample.wi;
                float p = sample.pdf;

                /* NEE */
                if(!sample.delta) {
                    LightSample light_sample = scene.light->connect_eye_path(sampler, intersection);
                    Vec3 s_dir = (light_sample.position-intersection.hitpoint);
                    float len = length(s_dir);
                    s_dir = s_dir/len;

                    Ray shadow_ray(intersection.hitpoint + s_dir*EPS,
                            s_dir,
                            len-2.0*EPS);


                    if(bsdf->non_zero(intersection,ray.d*(-1.0),shadow_ray.d)){
                    if(!scene.primitive->intersect_any(shadow_ray)){
                        float r = len;
                        col = col + mul*bsdf->eval(
                            ray.d*(-1.0),
                            shadow_ray.d,
                            intersection
                        )*light_sample.intensity/(light_sample.pdf*r*r);
                    }
                    }
                } else {
                    sampled_delta = true;
                }

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

                prev_intersection = intersection;
            }
        }
        return col;
    }

    float rF(const Vec3& radiance) const {
        return radiance.x+radiance.y+radiance.z;
    }

    float get_cdf(int i, const std::vector<float>& cdf) const {
        if(i<0){
            return 0.0;
        } else return cdf.at(i);
    }

    int choose_bootstrap(float ecs, const std::vector<float>& cdf) const {
        ecs = ecs*cdf.back() - 0.00001;
        //binary search
        int first = 0;
        int last = cdf.size();
        int middle;
        while(true){
            if((last-first)<5){
                for(int i = first; i<last-1; i++){
                    if(get_cdf(i-1,cdf)<=ecs && cdf.at(i)>=ecs){
                        return i;
                    }
                }
                return last-1;
            }

            middle = (last+first)/2;
            if(ecs<cdf.at(middle)){
                last = middle+1;
            } else {
                first = middle+1;
            }

        }
    }

    void add_image(RgbImage* im, float u, float v, Vec3 toAdd) const {
        int i = u*im->res.x;
        int j = (1.0-v)*im->res.y;

        Vec3* px = im->get_pixel(i,j);
                        //*px = Vec3(1.0,0.0,1.0);

        #pragma omp atomic
        (*px).x += toAdd.x;
        #pragma omp atomic
        (*px).y += toAdd.y;
        #pragma omp atomic
        (*px).z += toAdd.z;
    }

public:
    PathtracerMLT(const RenderSettings& rs): samples(rs.spp) {}

    void render(Scene& sc, std::string filename) {
        RgbImage* image= &(sc.camera->get_image());
        Vec2i res = image->res;

        int mChains = 32;
        int nBootstrap = res.x*res.y*8;
        std::vector<float> cdf(nBootstrap);
        std::vector<Vec3> values(nBootstrap);
        std::vector<std::vector<float>> record(nBootstrap);
        int nMutations = samples*samples*res.x*res.y/mChains;

        #pragma omp parallel for
        for(int i = 0; i<nBootstrap; i++){
            MLTSampler s(i);
            float r1 = s.sample();
            float r2 = s.sample();
            Vec3 h = integrate(sc, Vec2(r1,r2)*2.0-Vec2(1.0), s);
            cdf[i] = rF(h);
            values[i] = h;
            record.at(i) = s.sample_vector();
        }

        for(int i = 1; i<nBootstrap;i++){
            cdf[i] = cdf[i-1]+cdf[i];
        }

        float weight = float(res.x*res.y)/float(mChains*nMutations);
        float lcdf = cdf.back()/nBootstrap;



        #pragma omp parallel for
        for(int i = 0; i<mChains; i++){
            RandomSampler s(i);
            if(i%50==0)
                std::cout<<"rendering line "<<i<<std::endl;
        //for each chaic
            float ecs = s.sample();
            int j = choose_bootstrap(ecs,cdf);


            MLTSampler mlts(res.x+i*res.y +j);
            mlts.set_sample(record[j]);
            Vec3 value = values[j];
            float v = rF(value);
            add_image(image, mlts.sample_vector()[0],mlts.sample_vector()[1], (value/v)*lcdf*weight);


            for(int k = 1; k<nMutations; k++){
               // mlts.mutate();
               // add_image(image, mlts.old_sample_vector()[0],mlts.old_sample_vector()[1], (value/v)*lcdf*weight);
              //  mlts.reject();
              //  continue;
                mlts.mutate();
                float r1 = mlts.sample()*2.0-1.0;
                float r2 = mlts.sample()*2.0-1.0;
                Vec3 new_value = integrate(sc,Vec2(r1,r2),mlts);
                float new_v = rF(new_value);
                float alpha = new_v/v;
                //alpha = 0.0;
               // alpha = 1.1;
                if(new_v == 0.0){
                    alpha = 0.0;
                    new_v = 1.0;
                }
                if(alpha>1.0){
                    v = new_v;
                    value = new_value;
                    mlts.accept();
                    add_image(image, mlts.sample_vector()[0],mlts.sample_vector()[1], (value*(lcdf/v))*weight);
                } else {

                    add_image(image, mlts.old_sample_vector()[0],mlts.old_sample_vector()[1], (value*(lcdf/v))*(1.0-alpha)*weight);
                    add_image(image, mlts.sample_vector()[0],mlts.sample_vector()[1], (new_value*(lcdf)/new_v)*alpha*weight);
                    float r1 = s.sample();
                    if(r1<alpha){
                        v = new_v;
                        value = new_value;
                        mlts.accept();
                    } else {
                        mlts.reject();
                    }
                }
            }

        }


        image->save(filename);
    }
    static const char* name(){
        return "PathttracerMLT";
    }
};
}
#endif