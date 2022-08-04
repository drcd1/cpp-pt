#ifndef CPPPT_PATHTRACER
#define CPPPT_PATHTRACER

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
class Pathtracer : public Renderer{
    int samples;

    Vec3 render_sky(const Ray& r) const {
        return Vec3(0.0);
    }

    static float russian_roulette(const Vec3& col){
        //Todo: better rr
        return (col.x*0.2 + col.y*0.5 +col.z*0.3)*0.5 + 0.4;
    }



    Vec3 integrate(const Scene& scene, const Vec2& coords, Sampler& sampler) const {
        Ray ray = scene.camera->get_ray(coords);
        Intersection intersection;
        Vec3 col(0.0);
        Vec3 mul(1.0);
        for(int i = 0; i<32; i++){            
            bool intersected = scene.primitive->intersect(ray,&intersection);
            if(!intersected){
                col = col + mul*render_sky(ray);
                break;
            } else {
                col =  col + mul* intersection.material->emit(ray.d*(-1.0),intersection);
                Vec3 sample_direction;
                
                float p = intersection.material->sample(sampler, ray.d*(-1.0), intersection, &sample_direction);
                if(p<=0.0){
                    break;
                }

                Vec3 eval = intersection.material->eval(ray.d*(-1.0), sample_direction, intersection);
                if(i>2){
                    float rr = russian_roulette(eval);
                    if(sampler.sample() > rr) {
                        break;
                    }
                    p = p*rr;
                }

                mul = mul*eval/p;
                ray = Ray(intersection.hitpoint, sample_direction);
            }
        }
        return col;
    }

public:
    Pathtracer(int samples): samples(samples) {}

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
                acc = Vec3(linear2srgb(acc.x),linear2srgb(acc.y),linear2srgb(acc.z));
                image->put_pixel(i,j,acc);

            }
        }
        

        image->save(filename);
    }
};
}
#endif