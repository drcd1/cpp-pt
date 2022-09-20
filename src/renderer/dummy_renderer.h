#ifndef CPPPT_DUMMY_RENDERER
#define CPPPT_DUMMY_RENDERER

#include <renderer/renderer.h>
#include <math/math.h>
#include <camera/camera.h>
#include <image/rgb_image.h>
#include <primitive/ray.h>
#include <shape/intersection.h>
#include <scene.h>
#include <iostream>
namespace cpppt{
class DummyRenderer{
public:
    void render(Scene& sc, std::string filename) const {
        RgbImage* image= &(sc.camera->get_image());
        Vec2i res = image->res;
        for(int i = 0; i<res.x; i++){
            std::cout<<"rendering line "<<i<<std::endl;
            for(int j = 0; j<res.y; j++){
                Vector2<float> coords( ((float(i)+0.5)/float(res.x))*2.0-1.0, -(((float(j)+0.5)/float(res.y))*2.0-1.0));

                Ray ray = sc.camera->get_ray(coords);
                Intersection intersection;
                bool intersected = sc.primitive->intersect(ray,&intersection);

                #ifdef RAY_STATISTICS
                std::cout<<ray.tests<<std::endl;
                #endif
                if(!intersected){
                    image->put_pixel(i,j,Vec3(coords.x,coords.y,0.0));
                } else {
                    //float val = abs(dot(intersection.normal,ray.d));
                    //image->put_pixel(i,j,Vec3(val,val,val));
                    image->put_pixel(i,j, intersection.get_bxdf()->eval(ray.d*(-1.0),ray.d*(-1.0),intersection));
                }

                //std::cout<<coords.x<<" "<<coords.y<<std::endl;
                //


            }
        }

        image->save(filename);
    }
};
}
#endif