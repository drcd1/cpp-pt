#ifndef CPPPT_LIGHTTRACER
#define CPPPT_LIGHTTRACER

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
class Lighttracer : public Renderer{
    int samples;

    static float russian_roulette(const Vec3& col){
        //Todo: better rr
        return (fabs(col.x*0.2 + col.y*0.5 +col.z*0.3))*0.5 + 0.4;
    }


    void integrate(const Scene& scene, Sampler& sampler, float factor) const {



        LightPathStart lps = scene.light->sample(sampler);
        Ray ray(lps.position,lps.direction);

        Intersection intersection;
        Intersection prev_intersection;
        Vec3 col = lps.radiance;
        Vec3 mul = lps.radiance;
        bool sampled_delta = false;

        RgbImage& image = scene.camera->get_image();
        for(int i = 0; i<32; i++){
            bool intersected = scene.primitive->intersect(ray,&intersection);
            if(!intersected){

                break;
            } else {
                auto bsdf = intersection.material;

                Vec3 sample_direction;

                float p = bsdf->sample(sampler, ray.d*(-1.0), intersection, &sample_direction);
                Vec3 eval = bsdf->eval(ray.d*(-1.0), sample_direction, intersection);


                /* Connect to camera */
                if(!bsdf->is_delta()){
                    CameraConnection cc = scene.camera->connect_light_path(sampler, intersection);

                    Ray shadow_ray(intersection.hitpoint,
                        normalized(cc.pos-intersection.hitpoint),
                            length(cc.pos-intersection.hitpoint));

                    if(dot(shadow_ray.d,intersection.normal)>0.0 && cc.i >-1){
                    if(!scene.primitive->intersect_any(shadow_ray)){
                        Vec3 color =  mul*bsdf->eval(
                            ray.d*(-1.0),
                            shadow_ray.d,
                            intersection
                        ) * factor * 1.0/M_PI;

                        //color = color/ dot(shadow_ray.d, intersection.normal);

                        Vec3* px = image.get_pixel(cc.i,cc.j);
                        //*px = Vec3(1.0,0.0,1.0);



                        #pragma omp atomic
                        (*px).x += color.x;

                        #pragma omp atomic
                        (*px).y += color.y;
                        #pragma omp atomic
                        (*px).z += color.z;


                    }
                    }
                } else {
                    sampled_delta = true;
                }




                if(i>2){
                    float rr = russian_roulette(eval);
                    if(sampler.sample() > rr) {
                        break;
                    }
                    p = p*rr;
                }

                mul = mul*eval/p;
                ray = Ray(intersection.hitpoint, sample_direction);
                prev_intersection = intersection;
            }
        }


    }

public:
    Lighttracer(int samples): samples(samples) {}

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
                        integrate(sc,s,1.0/(samples*samples));
                    }
                }

            }
        }

        #pragma omp parallel for
        for(int i = 0; i<res.x; i++){
            for(int j = 0; j<res.y; j++){
                Vec3 h = *(image->get_pixel(i,j));
                h = Vec3(linear2srgb(h.x),linear2srgb(h.y),linear2srgb(h.z));
                image->put_pixel(i,j,h);
            }
        }



        image->save(filename);
    }
};
}
#endif