#ifndef CPPPT_LIGHTTRACER_MLT
#define CPPPT_LIGHTTRACER_MLT

#include <renderer/renderer.h>
#include <math/math.h>
#include <camera/camera.h>
#include <image/rgb_image.h>
#include <primitive/ray.h>
#include <shape/intersection.h>
#include <scene.h>
#include <math/sampler.h>
#include <bxdf/glossy_bxdf.h>
#include <omp.h>
#include <memory.h>
#include <iostream>
namespace cpppt{
class LighttracerMLT: public Renderer{
    int samples;


    struct RadianceSample{
        int i;
        int j;
        Vec3 radiance;
        RadianceSample(int i, int j, const Vec3& radiance):
        i(i),j(j),radiance(radiance){}
    };

    typedef std::vector<RadianceSample> RadianceSamples;

    static float russian_roulette(const Vec3& col){
        //Todo: better rr
        return std::min(1.0,(fabs(col.x*0.2 + col.y*0.5 +col.z*0.3))*0.5 + 0.4);

    }


    float rF(RadianceSamples& radiances) const {
        Vec3 radiance(0.0);
        for(auto v: radiances){
            radiance = radiance + v.radiance;
        }
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

    int max_path_length = 32;
    RadianceSamples integrate(const Scene& scene, Sampler& sampler) const {

        RadianceSamples ret;

        LightPathStart lps = scene.light->sample(sampler);
        Ray ray(lps.position,lps.direction);
        ray.o = ray.o+ray.d*EPS;

        Intersection intersection;
        Intersection prev_intersection;

        Vec3 mul = lps.radiance/lps.pdf;
        bool sampled_delta = false;

        RgbImage& image = scene.camera->get_image();
        for(int i = 0; i<max_path_length-1; i++){
            bool intersected = scene.primitive->intersect(ray,&intersection);
            if(!intersected){
                break;
            } else {
                auto bsdf = intersection.get_bxdf();
                if(bsdf->is_emitter()){
                    break;
                }

                DirectionalSample sample = bsdf->sample(sampler, ray.d*(-1.0), intersection);

                Vec3 sample_direction = sample.wi;
                float p = sample.pdf;

                Vec3 eval = bsdf->eval(ray.d*(-1.0), sample_direction, intersection);

                float correcting_factor = fabs(dot(intersection.normal,ray.d))/fabs(dot(intersection.g_normal,ray.d));


                /* Connect to camera */
                if(!sample.delta){
                    CameraConnection cc = scene.camera->connect_light_path(sampler, intersection);
                    Ray shadow_ray(intersection.hitpoint,
                        normalized(cc.pos-intersection.hitpoint),
                            length(cc.pos-intersection.hitpoint));
                    shadow_ray.o=shadow_ray.o+shadow_ray.d*EPS;

                    if(bsdf->non_zero(intersection,ray.d*(-1.0),shadow_ray.d) && cc.i >-1){
                    if(!scene.primitive->intersect_any(shadow_ray)){


                        Vec3 color =  mul*bsdf->eval(
                            normalized(ray.d*(-1.0)),
                            normalized(shadow_ray.d),
                            intersection
                        );// * factor;// * correcting_factor * fabs(dot(intersection.g_normal,shadow_ray.d)) / fabs(dot(intersection.normal,shadow_ray.d));

                        color = color*cc.factor;//(0.120/0.154);

                        Vec3* px = image.get_pixel(cc.i,cc.j);
                        //*px = Vec3(1.0,0.0,1.0);

                        ret.push_back({cc.i,cc.j,color});

                    }
                    }
                } else {
                    sampled_delta = true;
                }




                if(i>2){
                    float rr = russian_roulette(eval/p);
                    if(sampler.sample() > rr) {
                        break;
                    }
                    p = p*rr;
                }

                mul = mul*eval/p;//* correcting_factor * fabs(dot(intersection.g_normal,sample.wi)) / fabs(dot(intersection.normal,sample.wi));
                ray = Ray(intersection.hitpoint, sample_direction);
                ray.o=ray.o+ray.d*EPS;
                prev_intersection = intersection;
            }
        }
        return ret;
    }

    void add_image(RgbImage* im, int i, int j, Vec3 toAdd) const {
        Vec3* px = im->get_pixel(i,j);
                        //*px = Vec3(1.0,0.0,1.0);
        if(std::isnan(toAdd.x)||std::isnan(toAdd.y)||std::isnan(toAdd.z)){
            throw std::runtime_error("MORE NANS?!!");
        }
        #pragma omp atomic
        (*px).x += toAdd.x;
        #pragma omp atomic
        (*px).y += toAdd.y;
        #pragma omp atomic
        (*px).z += toAdd.z;
    }

public:
    LighttracerMLT(const RenderSettings& rs): samples(rs.spp) {}

    void render(Scene& sc, std::string filename) {
        RgbImage* image= &(sc.camera->get_image());
        Vec2i res = image->res;


        int mChains = 32;
        int nBootstrap = res.x*res.y;
        std::vector<float> cdf(nBootstrap);
        std::vector<RadianceSamples> values(nBootstrap);
        std::vector<std::vector<float>> record(nBootstrap);
        int nMutations = samples*samples*res.x*res.y/32;

        #pragma omp parallel for
        for(int i = 0; i<nBootstrap; i++){
            MLTSampler s(i);

            values[i] = integrate(sc, s);
            cdf[i] = rF(values[i]);
            record.at(i) = s.sample_vector();
        }

        for(int i = 1; i<nBootstrap;i++){
            cdf[i] = cdf[i-1]+cdf[i];
        }

        float lcdf = cdf.back()/nBootstrap;
        std::cout<<"lcdf: "<<lcdf<<std::endl;
       // lcdf = 1.0f;

        float weight = (float(res.x*res.y)/float(mChains*nMutations));

        #pragma omp parallel for
        for(int i = 0; i<mChains; i++){
            RandomSampler s(i);
            if(i%50==0)
                std::cout<<"rendering line "<<i<<std::endl;
            //for each chain

            float ecs = s.sample();
            int j = choose_bootstrap(ecs,cdf);
            MLTSampler mlts(res.x+i*res.y +j);
            mlts.set_sample(record[j]);
            auto value = values[j];
            float v = rF(value);
            if(v <= 1e-8)
                v=1e-8;




            for(auto rs: value){
                add_image(image, rs.i,rs.j,rs.radiance*lcdf/v*weight);
            }


            for(int k = 1; k<nMutations; k++){
                mlts.mutate();

                float r1 = mlts.sample()*2.0-1.0;
                float r2 = mlts.sample()*2.0-1.0;
                auto new_value = integrate(sc,mlts);
                float new_v = rF(new_value);
                float alpha = new_v/v;

                if(std::isnan(lcdf/v)||std::isnan(lcdf/new_v)){
                    throw(std::runtime_error("We have nans?!"));
                }

                float f1 = lcdf/v;
                float f2 = lcdf/new_v;
                if(v<1e-7){
                    f1 = 0.0;
                }
                if(new_v <1e-7){
                    alpha = 0.0;
                    f2 = 0.0;
                    new_v = 1e-7;
                }

                if(alpha>1.0){
                    v = new_v;
                    value = new_value;
                    mlts.accept();
                    for(auto rs: value){
                        add_image(image, rs.i,rs.j,rs.radiance*(f1)*weight);
                    }
                } else {
                    for(auto rs: value){
                        add_image(image, rs.i,rs.j,rs.radiance*(f1)*weight*(1.0-alpha));
                    }
                    for(auto rs: new_value){
                        add_image(image, rs.i,rs.j,rs.radiance*(f2)*weight*alpha);
                    }
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


        /*

        #pragma omp parallel for
        for(int i = 0; i<res.x; i++){
            for(int j = 0; j<res.y; j++){
                Vec3 h = *(image->get_pixel(i,j));
                h = Vec3(linear2srgb(h.x),linear2srgb(h.y),linear2srgb(h.z));
                image->put_pixel(i,j,h);
            }
        }
        */


        image->save(filename);
    }
    static const char* name(){
        return "LighttracerMLT";
    }
};
}
#endif