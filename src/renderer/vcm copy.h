#ifndef CPPPT_VCM
#define CPPPT_VCM



#include <renderer/renderer.h>
#include <math/math.h>
#include <camera/camera.h>
#include <image/rgb_image.h>
#include <primitive/ray.h>
#include <shape/intersection.h>
#include <primitive/primitive_leaf.h>
#include <scene.h>
#include <algorithms/kdtree.h>
#include <math/sampler.h>
#include <omp.h>
#include <iostream>
#include <light/shape_light.h>


#define MAX_PATH_LENGTH 32

namespace cpppt{
class VCM : public Renderer{

    pdfSum(const LightPath& from_camera, const LightPath& from_light){
        //Note: paths longer than MAX_PATH_LENGTH cannot connect
        

        //first step: fill in probability for forwardPDF for from_light and reversePDF for from_camera
        LightPath a = from_camera;
        LightPath b = from_light;

        LightPathVertex& lc = a.back();
        LightPathVertex& ll = b.back();


        Vec3 dir = lc.pt-ll.pt;
        float r = length(dir);
        dir=dir/r;

        //TODO: TRANSPOSED
        {
        
        float p = ll.bxdf->pdf(normalized(ll.previous->pt - ll->pt),dir,ll.it);
        lc.p_lt_acc = ll.p_lt_acc*fabs(dot(dir,lc.it.g_normal))/(r*r);

        }
        //NOT TRANSPOSED
        {        
        float p = ll.bxdf->pdf(normalized(lc.previous->pt - lc->pt),-dir,lc.it);
        lc.p_lt_acc = lc.p_pt_acc*fabs(dot(dir,ll.it.g_normal))/(r*r);
        }

        //Now all the paths probabilities are computed, and we can just sum


        //todo: okay, so this is way to many copies
        lc.insert(lc.end(),ll.rbegin(),ll.rend());

        return computePdfSumFullPath(lc);
    }
    
    //in the camera direction
    float computePdfSumFullPath(LightPath lc){
        float sum = 0.0;

        sum+=lc.at(0).p_lt_acc;
        for(int i = 1;i<lc.length()-1; i++){
            sum+=lc.at(i).p_pt_acc*lc.at(i+1).p_lt_acc;
        }
        sum+=lc.at(lc.length-1).pt_pt_acc;
    }


    struct Photon{
        Vec3 pt;
        //todo: remove repeated data
        Vec3 normal;
        Vec3 color;
        Vec3 dir;
    };

    struct LightPathVertex{
        Intersection it;
        std::shared_ptr<BxDF> bxdf;
        float p_lt_acc;
        float p_pt_acc;
        //note: p_pt might be offset;
        LightPathVertex* previous;
    };

    typedef std::vector<LightPathVertex> LightPath;

    std::vector<LightPath> lightpaths;



    KDTree<Photon> kdtree;

    float kernel_function(float d_sqr,float r_sqr) const {
        if(d_sqr>r_sqr)
            return 0.0;
        return 1.0/(M_PI*r_sqr);
    }

    int samples;

    Vec3 render_sky(const Ray& r) const {
        return Vec3(.0);
    }

    static float russian_roulette(const Vec3& col){
        //Todo: better rr
        return (fabs(col.x*0.2 + col.y*0.5 +col.z*0.3))*0.5 + 0.4;
    }


    LightPath emit_light_path(const Scene& scene, Sampler& sampler, float factor) {
        LightPath path;
        LightPathStart lps = scene.light->sample(sampler);
        Ray ray(lps.position,lps.direction);
        LightPathVertex lpv;
        //TODO
        lpv.g_normal = Vec3(0.0);
        lpv.p_lt_acc = lps.area_pdf;
        lpv.pt = lps.position;
        lpv.p_pt_acc = 1.0;
        lpv.previous = nullptr;
        lpv.bxdf = nullptr; //todo
        //lpv.intersection = nullptr;//todo

        path.push_back(lpv);


        Intersection intersection;
        Intersection prev_intersection;
        Vec3 mul = lps.radiance/lps.pdf;

        RgbImage& image = scene.camera->get_image();
        float old_p = 1.0;
        for(int i = 0; i<MAX_PATH_LENGTH; i++){
            bool intersected = scene.primitive->intersect(ray,&intersection);
            if(!intersected){

                break;
            } else {
                auto bsdf = intersection.get_bxdf();
                if(bsdf->is_emitter()){
                    break;
                }

                LightPathVertex lpv;
                lpv.it = intersection;
                lpv.p_lt_acc = 1.0;
                if(i == 0){
                    lpv.p_pt_acc = scene.camera->pdf(intersection.hitpoint,intersection.g_normal);
                } else {
                    lpv.p_pt_acc = path.back().p_pt_acc*old_p;
                }
                path.push_back(lpv);


                DirectionalSample sample = bsdf->sample(sampler, ray.d*(-1.0), intersection);

                Vec3 sample_direction = sample.wi;
                float p = sample.pdf;

                Vec3 eval = bsdf->eval(ray.d*(-1.0), sample_direction, intersection);

                float correcting_factor = fabs(dot(intersection.normal,ray.d))/fabs(dot(intersection.g_normal,ray.d));

                // Store photon
                if(!sample.delta){
                    Photon p;
                    p.pt = intersection.hitpoint;
                    p.normal = intersection.normal;
                    p.color = mul*factor;
                    p.dir = ray.d*(-1.0);
                    //TODO: correcting factor

                    #pragma omp critical
                    {
                    kdtree.add(p);
                    }
                }

                if(!sample.delta){
                    CameraConnection cc = scene.camera->connect_light_path(sampler, intersection);

                    Ray shadow_ray(intersection.hitpoint,
                        normalized(cc.pos-intersection.hitpoint),
                            length(cc.pos-intersection.hitpoint));

                    if(bsdf->non_zero(intersection,ray.d*(-1.0),shadow_ray.d) && cc.i >-1){
                    if(!scene.primitive->intersect_any(shadow_ray)){


                        Vec3 color =  mul*bsdf->eval(
                            normalized(ray.d*(-1.0)),
                            normalized(shadow_ray.d),
                            intersection
                        ) * factor * correcting_factor * fabs(dot(intersection.g_normal,shadow_ray.d)) / fabs(dot(intersection.normal,shadow_ray.d));

                        color = color*cc.factor;

                        //note: cameras have lens with zero area
                        float p = path.back().p_lt_acc;
                        float p_forward = scene.camera->pdf(path.back().pt, path.back().g_normal);
                        float denom = p;
                        float connection_prob = factor; // 1/N_LIGHT_PATHS
                        //Todo: replace
                        for(int k = path.size()-1; k>0; k--){
                            //when computing the weights, k is the last vertex generated by pathtracing
                            //k = 0 is pathtracing, not counted in this loop
                            //k = 1 is nee
                            //k = path.size() is lighttracing, not counted in this loop
                            p_forward*=path.at(k).p_pt_acc;
                            denom += p_forward*path.at(k).p_pt_acc*connection_prob*path.at(k).p_lt_acc;
                        }
                        denom += p_forward*path.at(0).p_pt_acc;
                        float w = p/denom;
                        color = color*w;


                        Vec3* px = image.get_pixel(cc.i,cc.j);

                        #pragma omp atomic
                        (*px).x += color.x;

                        #pragma omp atomic
                        (*px).y += color.y;

                        #pragma omp atomic
                        (*px).z += color.z;


                    }
                    }
                }


                auto old_v = path.back();
                float r2 = lensqr(old_v.pt - intersection.hitpoint);
                float c = fabs(dot(old_v.g_normal,ray.d));
                //TODO: russian roulette prob
                old_v.p_pt_acc*= bsdf->pdf(sample.wi, ray.d*(-1.0),intersection)*r2/c;


                if(i>2){
                    float rr = russian_roulette(eval);
                    if(sampler.sample() > rr) {
                        break;
                    }
                    p = p*rr;
                }
                old_p = p;

                mul = mul*eval/p* correcting_factor * fabs(dot(intersection.g_normal,sample.wi)) / fabs(dot(intersection.normal,sample.wi));

                ray = Ray(intersection.hitpoint, sample_direction);
                prev_intersection = intersection;
            }
        }


    }

    Vec3 integrate(const Scene& scene, const Vec2& coords, Sampler& sampler) const {
        Ray ray = scene.camera->get_ray(coords);
        Intersection intersection;
        Vec3 col(0.0);
        Vec3 mul(1.0);
        bool sampled_delta = true;
        float p = 1.0;
        LightPath path;


        std::vector<LightPathVertex> lc;
        LightPathVertex lv;
        //TODO: init camera vertex
        lc.push_back(lv);

        for(int i = 0; i<MAX_PATH_LENGTH; i++){
            bool intersected = scene.primitive->intersect(ray,&intersection);



            if(!intersected){
                //col = col + mul*render_sky(ray)*compute_pt_direct_factor();
                break;
            } else {
                auto bsdf = intersection.get_bxdf();
                LightPathVertex lv;
                lv.p_pt_acc = p*dot();
                lc.push_back()

                if(bsdf->is_emitter()){
                    if(sampled_delta){

                        //scene.light = dynamic_cast
                        //float p_bsdf_nee = scene.light->pdf(static_cast<std::shared_ptr<EmissiveBxDF>>(bsdf))
                        //float p_bsdf_nee  = scene.light->pdf()*(r*r)/bsdf_cos;

                        Vec3 rad = bsdf->emit(ray.d*(-1.0),intersection);

                        col = col + mul*rad;
                    } else {
                        

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


                /* BDPTConnect */
                if(!sample.delta){
                    std::cout<<"delta"<<std::endl;
                   /* for(lv: p){
                        light_sample = lv;

                    }*//*
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
                    }*/
                    //}
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
                //mul = mul*eval/p;
                //ray = Ray(intersection.hitpoint+sample_direction*EPS, sample_direction);
            }
        }
        return col;
    }

public:
    VCM(const RenderSettings& rs): samples(rs.spp) {}

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
                image->put_pixel(i,j,acc);

            }
        }
        }


        image->save(filename);
    }
};
}
#endif