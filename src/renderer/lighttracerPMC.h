#ifndef CPPPT_LIGHTTRACER_PMC
#define CPPPT_LIGHTTRACER_PMC

#include <renderer/renderer.h>
#include <math/math.h>
#include <camera/camera.h>
#include <image/rgb_image.h>
#include <primitive/ray.h>
#include <shape/intersection.h>
#include <scene.h>
#include <math/sampler.h>
#include <light/light.h>
#include <bxdf/glossy_bxdf.h>
#include <omp.h>
#include <memory.h>
#include <iostream>
namespace cpppt{
class LighttracerPMC: public Renderer{
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
        return std::max(2.0f,radiance.x+radiance.y+radiance.z);
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
      
                    if(bsdf->name != "GlossyBxDF" && bsdf->non_zero(intersection,ray.d*(-1.0),shadow_ray.d) && cc.i >-1){
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
            std::cout<<"NANS MORe!?"<<std::endl;
            return;
            
            //throw std::runtime_error("MORE NANS?!!");
        }
        #pragma omp atomic
        (*px).x += toAdd.x;
        #pragma omp atomic
        (*px).y += toAdd.y;
        #pragma omp atomic
        (*px).z += toAdd.z;
    }


    std::vector<RandomSampler> uniform_samplers;

public:
    LighttracerPMC(const RenderSettings& rs): samples(rs.spp) {}

    void render(Scene& sc, std::string filename) {
        RgbImage* image= &(sc.camera->get_image());
        Vec2i res = image->res;



        struct PopulationSample{
            PMCSampler s;
            //a value proportional to the pdf
            float v;
            //a value proportional to the pdf, taking into account the resampling process
            float target;
            RadianceSamples value;
        };

        //TODO: check this weight

        int popSize = res.x*res.y;
        int nIter = samples*samples;
        float weight = (float(res.x*res.y)/float(popSize*nIter));


        //note: this is a mutation over the population
        int nMutations = samples*samples*res.x*res.y/popSize;


        std::vector<PopulationSample> population(popSize);

        for(int i = 0; i<omp_get_max_threads(); i++){
            uniform_samplers.push_back(RandomSampler(i));
            //evaluate the samples
            //add to the image
        }
        
        #pragma omp parallel for
        for(int i = 0; i<popSize; i++){
            int tid = omp_get_thread_num();
            population.at(i).s = PMCSampler(&uniform_samplers[tid]);
            population.at(i).value = integrate(sc, population.at(i).s);
            
        }



        
        float max = 0.0f;
        float targetSum = 0.0f;
        float avg = 1.0f;
        for(int i = 0; i<popSize; i++){
            float r = rF(population.at(i).value);
            max = r>max?r:max;
            targetSum+=r;
            population.at(i).target = r;
            population.at(i).v = r;
            population.at(i).value=population.at(i).value;
            //for(auto rs: population.at(i).value){
                //todo: does this makes sense?
                //add_image(image, rs.i,rs.j,rs.radiance*weight);
            //}
        }
        //targetSum = 1.0f;
        bool resample = true;
        bool markov = false;

        avg=targetSum*avg/popSize;
        RandomSampler s(1024);
        for(int i = 0; i<nMutations; i++){
            if(resample || i==1){
                //avg=targetSum*avg/popSize;
                std::vector<PopulationSample> new_pop;
                std::copy_if(population.begin(), population.end(), std::back_inserter(new_pop),[&](PopulationSample& i){
                    return s.sample()<i.target/max;
                });

                int size = new_pop.size();
                

                new_pop.reserve(popSize);
                for(int i = new_pop.size(); i<popSize; i++){
                    PopulationSample ps;
                    
                    int j = int(s.sample()*size)%size;
                    ps.target=new_pop.at(j).target;
                    ps.s = new_pop.at(j).s;
                    ps.v = new_pop.at(j).v;
                    ps.value = new_pop.at(j).value;

                    new_pop.push_back(ps);

                }


                std::swap(population, new_pop);   
            }

            #pragma omp parallel for
            for(int i = 0; i<popSize;i++){
                population.at(i).s.set_sampler(&uniform_samplers.at(omp_get_thread_num()));
                population.at(i).s.mutate();

                auto new_value = integrate(sc,population.at(i).s);

                float new_v = rF(new_value);
                if(markov){
                if(new_v>=population.at(i).v){
                    float f1= avg/(population.at(i).v);
                
                    population.at(i).s.accept();
                    population.at(i).value = new_value;
                    population.at(i).v = new_v;
                    population.at(i).target = new_v/population.at(i).v;
                    for(auto rs: new_value){
                        //todo: does this makes sense?
                        add_image(image, rs.i,rs.j,rs.radiance*weight*f1);
                    }

                } else {
                    float alpha = new_v/population.at(i).v;

                    float f1= avg/(population.at(i).v);
                    float f2= avg/(new_v);

                    for(auto rs: population.at(i).value){
                        //todo: does this makes sense?
                        add_image(image, rs.i,rs.j,rs.radiance*weight*(1.0f-alpha)*f1);
                    }
                    for(auto rs: new_value){
                        //todo: does this makes sense?
                        add_image(image, rs.i,rs.j,rs.radiance*weight*(alpha)*f2);
                    }


                    //accept or reject
                    float ecs = uniform_samplers[omp_get_thread_num()].sample();
                    if(ecs<alpha){
                        //accept
                        population.at(i).target = new_v/population.at(i).v;
                        population.at(i).s.accept();
                        population.at(i).value = new_value;
                        population.at(i).v = new_v;
                    } else {
                        population.at(i).s.reject();
                        population.at(i).target = 1.0f;
                    }
                }
                } else {

                    //target = new_v/sum*old_sum/v

                    float f1 = avg/population.at(i).v;
                    for(auto rs: new_value){
                        add_image(image, rs.i,rs.j,rs.radiance*weight*f1);
                    }
                    population.at(i).target = new_v/population.at(i).v;
                    population.at(i).s.accept();
                    population.at(i).value = new_value;
                    population.at(i).v = new_v;

                }

            }



            if(resample){
                max = 0.0f;
                targetSum = 0.0f;
                for(int i = 0; i<popSize; i++){
                    float r = population.at(i).target;
                    max = r>max?r:max;
                    targetSum+=population.at(i).target;
                }
            }
            
            

        }
        

        image->save(filename);
    }
    static const char* name(){
        return "LighttracerPMC";
    }
};
}
#endif