#ifndef CPPPT_PPG_LT
#define CPPPT_PPG_LT

#include <renderer/renderer.h>
#include <math/math.h>
#include <camera/camera.h>
#include <image/rgb_image.h>
#include <primitive/ray.h>
#include <algorithms/sdtree.h>
#include <shape/intersection.h>
#include <scene.h>
#include <math/sampler.h>
#include <omp.h>
#include <iostream>

namespace cpppt{


class PPG_LT : public Renderer{
    int samples;

    SDTree* st;

    Vec3 render_sky(const Ray& r) const {
        return Vec3(.0);
    }

    static float russian_roulette(const Vec3& col){
        //Todo: better rr
        return (fabs(col.x*0.2 + col.y*0.5 +col.z*0.3))*0.5 + 0.4;
    }

    struct Vertex{
        Vec3 pos;
        Vec3 dir;
        Vec3 radiance;
    };
    float rF(const Vec3& radiance) const {
        return std::max(radiance.x,std::max(radiance.y,radiance.z));
    }
    Vec3 integrate(const Scene& scene, const Vec2& coords, Sampler& sampler) {
        Ray ray = scene.camera->get_ray(coords);
        Intersection intersection;
        Intersection prev_intersection;
        Vec3 col(0.0);
        Vec3 mul(1.0);
        std::vector<Vertex> vertices;

        for(int i = 0; i<32; i++){
            bool intersected = scene.primitive->intersect(ray,&intersection);
            if(i>0){
                vertices.push_back(Vertex());
                vertices.back().pos = ray.o;
                vertices.back().dir = ray.d;
                vertices.back().radiance = Vec3(1.0);
            }
            if(!intersected){
                Vec3 rad = render_sky(ray);
                col = col + mul*rad;
                if(i>0){
                    vertices.back().radiance = rad;
                }
                break;
            } else {
                auto bsdf = intersection.get_bxdf();
                if(bsdf->is_emitter()){
                    Vec3 rad = bsdf->emit(ray.d*(-1.0),intersection);
                    if(i>0){
                        vertices.back().radiance = rad;
                    }
                    col = col + mul*rad;
                    break;
                }
                DirectionalSample sample;
                float weight = 1.0;
                if(sampler.sample()>0.5){
                    sample = bsdf->sample(sampler, ray.d*(-1.0), intersection);
                    weight = (sample.pdf/(sample.pdf + st->pdf(intersection.hitpoint, sample.wi)));
                    if(std::isnan(weight) || sample.pdf==0){
                        std::cout<<"nan";
                    }
                } else {



                    sample = st->sample(sampler, intersection.hitpoint);
                    if(!bsdf->non_zero(intersection,ray.d*(-1.0),sample.wi)){
                        break;
                    }
                    float bsdf_pdf =  bsdf->pdf(ray.d*(-1.0), sample.wi, intersection);
                    weight = (sample.pdf/(sample.pdf +bsdf_pdf));
                    if(std::isnan(weight)|| sample.pdf == 0.0){
                        std::cout<<"nan";
                    }

                }

                if(std::isnan(sample.pdf) || std::isnan(weight) || std::isnan(weight/sample.pdf)){
                    std::cout<<"NaNi?"<<std::endl;
                }

                Vec3 sample_direction = sample.wi;
                float p = sample.pdf*0.5;

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

                //note: this is the local multiplication factor of the next vert.
                Vec3 new_mul =  eval*(weight/p);
                if(i>0){
                    vertices.back().radiance =new_mul;
                }
                mul = mul*new_mul;
                ray = Ray(intersection.hitpoint+sample_direction*EPS, sample_direction);
                prev_intersection = intersection;

            }
        }
        //only splat if radiance is enough

        if(vertices.size()>0 &&(vertices.back().radiance.x>0 || vertices.back().radiance.y>0 || vertices.back().radiance.z>0)){
            for(int i = vertices.size()-1; i>=0; i--){
                if(i<vertices.size()-1){
                    vertices.at(i).radiance = vertices.at(i).radiance*vertices.at(i+1).radiance;
                }

                auto& v = vertices.at(i);
                st->splat(rF(v.radiance),v.pos,v.dir);
            }
        }

        return col;
    }

public:
    PPG_LT(const RenderSettings& rs): samples(rs.spp), st(nullptr){}

    void render(Scene& sc, std::string filename) {
        st = new SDTree({sc.primitive->get_bounds().min,sc.primitive->get_bounds().max});
        RgbImage* image= &(sc.camera->get_image());
        Vec2i res = image->res;
        int total_samples = samples*samples;
        std::vector<int> sample_vector;
        int spp = 1;
        while(true){
            sample_vector.push_back(spp);
            spp*=2;
            if(spp>=total_samples){
                break;
            }
        }

        for(int l=0; l<sample_vector.size(); l++){
            std::ostringstream oss;
            oss<<l;
            st->debug_log("iteration_"+oss.str());
            if(l==sample_vector.size()-1){
                st->no_learning();
            }
            int spp = sample_vector.at(l);

            #pragma omp parallel for
            for(int i = 0; i<res.x; i++){
                RandomSampler s(i);
                if(i%50==0)
                    std::cout<<"rendering line "<<i<<std::endl;
                for(int j = 0; j<res.y; j++){
                    Vec3 acc(0.0);
                    for(int k = 0; k<spp; k++){
                        float r1 = s.sample();
                        float r2 = s.sample();

                        Vector2<float> coords( ((float(i)+r1)/float(res.x))*2.0-1.0, -(((float(j)+r2)/float(res.y))*2.0-1.0));
                        acc = acc +integrate(sc, coords, s);


                    }
                    acc = acc/(float(spp));
                    if(l==sample_vector.size()-1)
                        image->put_pixel(i,j,acc);

                }
            }
            /*
            if(l==sample_vector.size()-1){
                float est = 0.0;
                RandomSampler s(2);
                for(int i = 0; i<100000; i++){
                    auto sd = st->sample(s, Vec3(0.4,0.4,0.0));
                    float p = sd.pdf;
                    est = (est*i + 1.0/(4.0*M_PI*p))/float(i+1);
                    if(i%1000==0){
                        std::cout<<"Current Estimate: "<<est<<std::endl;
                    }
                }
            }*/
            st->update();
        }

        delete st;

        image->save(filename);
    }

    static const char* name(){
        return "PPG LT";
    }
};
}
#endif