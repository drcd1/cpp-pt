#ifndef CPPPT_DEBUG_RENDERER
#define CPPPT_DEBUG_RENDERER

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
class DebugRenderer : public Renderer{

    Vec3 render_sky(const Ray& r) const {
        return Vec3(1.0,0.0,1.0);
    }


    Vec3 integrate(const Scene& scene, const Vec2& coords, Sampler& sampler) const {
        Ray ray = scene.camera->get_ray(coords);
        Intersection intersection;
        Vec3 col(0.0);
        Vec3 mul(1.0);
        bool intersected = scene.primitive->intersect(ray,&intersection);
        if(!intersected){
            return render_sky(ray);
        } else {
            return Vec3(abs(dot(intersection.normal, -ray.d)));
        }

    }

public:

    void render(Scene& sc, std::string filename) {
        RgbImage* image= &(sc.camera->get_image());
        Vec2i res = image->res;

        #pragma omp parallel for
        for(int i = 0; i<res.x; i++){
            RandomSampler s(i);
            if(i%50==0)
                std::cout<<"rendering line "<<i<<std::endl;
            for(int j = 0; j<res.y; j++){
                Vec3 acc(0.0);
                Vector2<float> coords( ((float(i))/float(res.x))*2.0-1.0, -(((float(j))/float(res.y))*2.0-1.0));
                acc = acc +integrate(sc, coords, s);


                image->put_pixel(i,j,acc);

            }
        }


        image->save(filename);
    }
};
}
#endif