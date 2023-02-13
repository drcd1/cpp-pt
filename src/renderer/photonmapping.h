#ifndef CPPPT_PHOTONMAPPING
#define CPPPT_PHOTONMAPPING

#include <renderer/renderer.h>
#include <math/math.h>
#include <camera/camera.h>
#include <image/rgb_image.h>
#include <primitive/ray.h>
#include <algorithms/kdtree.h>
#include <shape/intersection.h>
#include <scene.h>
#include <math/sampler.h>
#include <omp.h>
#include <iostream>
namespace cpppt{
class PhotonMapping : public Renderer{
    int samples;

    static float russian_roulette(const Vec3& col){
        //Todo: better rr
        return (fabs(col.x*0.2 + col.y*0.5 +col.z*0.3))*0.5 + 0.4;
    }
    Vec3 render_sky(const Ray& r) const {
        return Vec3(.0);
    }

    struct Photon{
        Vec3 pt;
        //todo: remove repeated data
        Vec3 normal;
        Vec3 color;
        Vec3 dir;
    };

    KDTree<Photon> kdtree;

    float kernel_function(float d_sqr,float r_sqr) const {
        if(d_sqr>r_sqr)
            return 0.0;
        return 1.0/(M_PI*r_sqr);
    }

    Vec3 integrate(const Scene& scene, const Vec2& coords, Sampler& sampler) {
        Ray ray = scene.camera->get_ray(coords);
        Intersection intersection;
        Vec3 col(0.0);
        Vec3 mul(1.0);
        {
            bool intersected = scene.primitive->intersect(ray,&intersection);
            if(!intersected){
                col = col + mul*render_sky(ray);
                return col;
            } else {
                auto bsdf = intersection.get_bxdf();
                if(bsdf->is_emitter()){
                    col = col + mul*bsdf->emit(ray.d*(-1.0),intersection);
                    return col;
                }
                float radius = 0.01;
                float dist2 = lensqr(intersection.hitpoint - ray.o);

                kdtree.run_kernel(intersection.hitpoint,radius,[&](int a){
                    float w = kernel_function(lensqr(intersection.hitpoint- kdtree.points.at(a).pt),radius*radius);
                    w = dot(intersection.normal, kdtree.points.at(a).normal)>0.0?w:0.0;
                    if(w<0.001) return;
                    w = w/fabs(dot(kdtree.points.at(a).normal, kdtree.points.at(a).dir));
                    col=col+kdtree.points.at(a).color*bsdf->eval(ray.d*(-1.0), kdtree.points.at(a).dir, intersection)*w;//*0.001;
                });


            }
        }
        return col;
    }

    void emit_light_path(const Scene& scene, Sampler& sampler, float factor) {



        LightPathStart lps = scene.light->sample(sampler);
        Ray ray(lps.position,lps.direction);

        Intersection intersection;
        Intersection prev_intersection;
        Vec3 mul = lps.radiance/lps.pdf;

        RgbImage& image = scene.camera->get_image();
        for(int i = 0; i<32; i++){
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


                /* Store photon */
                if(!sample.delta){
                    Photon p;
                    p.pt = intersection.hitpoint;
                    p.normal = intersection.normal;
                    p.color = mul*factor;
                    p.dir = ray.d*(-1.0);

                    #pragma omp critical
                    {
                    kdtree.add(p);
                    }
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
    PhotonMapping(const RenderSettings& rs): samples(rs.spp) {}

    void render(Scene& sc, std::string filename) {
        RgbImage* image= &(sc.camera->get_image());
        Vec2i res = image->res;

        for(int k =0; k<samples*samples; k++){
        #pragma omp parallel for
        for(int i = 0; i<res.x; i++){



            RandomSampler s(k*res.x+i);
            for(int j = 0; j<res.y; j++){
                Vec3 acc(0.0);
                emit_light_path(sc,s,1.0/(res.x*res.y));

            }
        }

        kdtree.build();

        #pragma omp parallel for
        for(int i = 0; i<res.x; i++){


            RandomSampler s(k*res.x+i);
            for(int j = 0; j<res.y; j++){
                Vec3 acc(0.0);
                float r1 = s.sample();
                float r2 = s.sample();

                Vector2<float> coords( ((float(i)+r1)/float(res.x))*2.0-1.0, -(((float(j)+r2)/float(res.y))*2.0-1.0));
                acc = acc +integrate(sc, coords, s);



                acc = acc/(float(samples*samples));
                image->put_pixel(i,j,*(image->get_pixel(i,j))+acc);

            }
        }
        kdtree = KDTree<Photon>();

        std::cout<<"rendered sample "<<k<<std::endl;
        }

        image->save(filename);
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


    }

};
}
#endif