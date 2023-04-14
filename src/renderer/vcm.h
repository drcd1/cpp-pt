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


    struct Photon{
        Vec3 pt;
        //todo: remove repeated data
        Vec3 normal;
        Vec3 color;
        Vec3 dir;
    };

    struct LightPathVertex{


        Intersection it;
        union {
            const BxDF* bxdf = nullptr;
            const Light* light;
        };
        bool is_light = false;

        float p_lt_acc;
        float p_pt_acc;
        Vec3 radiance = Vec3(0.0f);

    };

    typedef std::vector<LightPathVertex> LightPath;

    std::vector<LightPath> lightpaths;



    KDTree<Photon> kdtree;
/*

    bool testDL(){
        //create simple path
        //create the path with: LT, PT and VC
        //check if probs computed are the same
        //weights are the same
        //they sum to one

        LightPath lp;
        LightPath cp;

        Intersection vc;
        Intersection vmid;
        Intersection vl;

        vc.hitpoint = Vec3(0.f,0.f,0.f);
        vmid.hitpoint = Vec3(0.0f,2.0f, -1.0f);
        vmid.g_normal = Vec3(0.0f,0.0f,1.0f);
        vl.hitpoint = Vec3(0.0f,0.0f,1.0f);
        vl.g_normal = Vec3(0.0f,-1.0f,0.0f);

        {

        LightPathVertex v1;
      //  v1.it = v1.hitpoint;
        LightPathVertex v2;
        LightPathVertex v3;

        lp.push_back(v1);
        lp.push_back(v2);
        lp.push_back(v3);

        }

        {

        LightPathVertex v1;
        LightPathVertex v2;
        LightPathVertex v3;

        lc.push_back(v1);
        lc.push_back(v2);
        lc.push_back(v3);

        }




    }
*/

    float compute_area_pdf(const LightPathVertex& v, const LightPathVertex& v_from, const Vec3& from_from, bool transposed) const {
        //TODO: account for transposed
        Vec3 wo = normalized(from_from-v_from.it.hitpoint);
        Vec3 wi = v.it.hitpoint-v_from.it.hitpoint;

        float r= length(wi);
        wi = wi/r;

        float cosTheta = fabs(dot(wi,v.it.g_normal));

        float p = v_from.bxdf->pdf(wo,wi,v_from.it);
        float rr = russian_roulette(v_from.bxdf->eval(wo,wi,v_from.it));

        return compute_area_pdf(p,r*r,cosTheta);
    }
    float compute_area_pdf(const LightPathVertex& v, const LightPathVertex& v_from, float solid_angle_pdf) const {
        //TODO: account for transposed
        Vec3 wi = v.it.hitpoint-v_from.it.hitpoint;

        float r= length(wi);
        wi = wi/r;
        float cosTheta = fabs(dot(wi,v.it.g_normal));

        return compute_area_pdf(solid_angle_pdf,r*r,cosTheta);
    }

    inline float compute_area_pdf_r(float solid_angle_pdf, float r, float cosTheta) const {
        return compute_area_pdf(solid_angle_pdf,r*r,cosTheta);
    }

    inline float compute_area_pdf(float solid_angle_pdf, float rSQ, float cosTheta) const {
        return solid_angle_pdf*cosTheta/rSQ;
    }


    //assumes every vertex has the conditional
    //probabilities expressed in terms of area
    void accPdf(LightPath& lp) const {
        for(int i = 1; i<lp.size(); i++){
            lp.at(i).p_pt_acc *= lp.at(i-1).p_pt_acc;
            lp.at(lp.size()-i-1).p_lt_acc *= lp.at(lp.size()-i).p_lt_acc;
        }
    }


    LightPath pdfSum(const LightPath& from_camera, const LightPath& from_light, float& pdfSum) const {
        //Note: paths longer than MAX_PATH_LENGTH cannot connect
        LightPath l;

        if(from_light.size()>0 && from_camera.size()>0){

        //first step: fill in probability for forwardPDF for from_light and reversePDF for from_camera
        LightPath a = from_camera;
        LightPath b = from_light;

        LightPathVertex& lc = a.back();
        LightPathVertex& ll = b.back();


        Vec3 dir = lc.it.hitpoint-ll.it.hitpoint;
        float r = length(dir);
        dir=dir/r;

        //TODO: TRANSPOSED
        if(b.size()>1)
        {
        Vec3 wo = normalized(b.at(b.size()-2).it.hitpoint - ll.it.hitpoint);
        Vec3 wi = dir;

        float p = ll.bxdf->pdf(wo,wi,ll.it);
        p*= russian_roulette(ll.bxdf->eval(wo,wi,ll.it));
        //lc.p_lt_acc = ll.p_lt_acc*fabs(dot(dir,lc.it.g_normal))/(r*r);
        float cos_theta = fabs(dot(dir,lc.it.g_normal));
        lc.p_lt_acc =  p*cos_theta/(r*r);

        } else {
            Vec3 wi = dir;
            //TODO: don't hardcode!
            float p = ll.p_lt_acc*fabs(dot(wi,ll.it.g_normal))/M_PI;
            float cos_theta = fabs(dot(wi,lc.it.g_normal));
            lc.p_lt_acc = p*cos_theta/(r*r);
        }

        //NOT TRANSPOSED
        {
        Vec3 wi = dir*(-1.0);
        Vec3 wo = normalized(a.at(a.size()-2).it.hitpoint - lc.it.hitpoint);

        float p = lc.bxdf->pdf(wo,wi,lc.it);
        float cos_theta = fabs(dot(dir,ll.it.g_normal));
        p*= russian_roulette(lc.bxdf->eval(wo,wi,lc.it));
        p*=cos_theta/(r*r);

        ll.p_pt_acc = p;
        }

        //Now all the paths probabilities are computed, and we can just sum


        //todo: okay, so this is way to0 many copies
        l = a;
        l.insert(l.end(),b.rbegin(),b.rend());
        } else if(from_camera.size()>0){
            l = from_camera;
        } else {
            l.insert(l.end(), from_light.rbegin(), from_light.rend());
        }
        accPdf(l);

        pdfSum = computePdfSumFullPath(l);
        return l;
    }

    //in the camera2light  direction
    float computePdfSumFullPath(const LightPath& lc) const {
        float sum = 0.0;

        sum+=lc.at(0).p_lt_acc;
        for(int i = lc.size()-2;i<lc.size()-1; i++){
            sum+=lc.at(i).p_pt_acc*lc.at(i+1).p_lt_acc;
            if(std::isnan(sum)){
                std::cout<<"ON NOO"<<std::endl;
            }
        }
        sum+=lc.at(lc.size()-1).p_pt_acc;
        if(std::isnan(sum)){
            std::cout<<"klkadlkw"<<std::endl;
        }
        return sum;
    }


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
        return std::min(1.0f,(fabs(col.x*0.2f + col.y*0.5f +col.z*0.3f))*0.5f + 0.4f);
    }


    LightPath emit_light_path(const Scene& scene, Sampler& sampler, float factor) {
        LightPath path;
        LightPathStart lps = scene.light->sample(sampler);
        Ray ray(lps.position,lps.direction);
        ray.o = ray.o+ray.d*EPS;
        LightPathVertex lpv;

        lpv.it = Intersection();
        lpv.it.hitpoint = lps.position;
        lpv.it.g_normal = lps.normal;
        lpv.p_lt_acc = lps.area_pdf;
        lpv.p_pt_acc = 1.0;
        //MISSING COSINE FACTOR COMPUTATION?! HOW TO DISTINGUISH DELTA FROM NOT?
        //irradiance    only for first vertex
        //should factor be 1.0?
        lpv.radiance = lps.radiance;
        lpv.light = lps.light; //todo
        lpv.is_light = true;
        //lpv.intersection = nullptr;//todo

        path.push_back(lpv);


        Intersection intersection;
        Intersection prev_intersection;
        Vec3 mul = lps.radiance/lps.pdf;

        RgbImage& image = scene.camera->get_image();
        float old_p = lps.angle_pdf;
        for(int i = 0; i<32; i++){
            bool intersected = scene.primitive->intersect(ray,&intersection);
            if(!intersected){
                return path;
            } else {
                auto bsdf = intersection.get_bxdf();
                if(bsdf->is_emitter()){
                    return path;
                }

                LightPathVertex lpv;
                lpv.it = intersection;
                lpv.bxdf = bsdf.get();

                {

                lpv.p_lt_acc = compute_area_pdf(lpv,path.back(),old_p);
                }

                lpv.radiance = mul;
                path.push_back(lpv);



                DirectionalSample sample = bsdf->sample(sampler, ray.d*(-1.0), intersection);

                Vec3 sample_direction = sample.wi;
                float p = sample.pdf;

                Vec3 eval = bsdf->eval(ray.d*(-1.0), sample_direction, intersection);

                float correcting_factor = fabs(dot(intersection.normal,ray.d))/fabs(dot(intersection.g_normal,ray.d));

                // Store photon
                /*
                if(!sample.delta){
                    Photon p;
                    p.pt = intersection.hitpoint;
                    p.normal = intersection.normal;
                    p.color = mul;
                    p.dir = ray.d*(-1.0);
                    //TODO: correcting factor

                    #pragma omp critical
                    {
                    kdtree.add(p);
                    }
                }
                */

                if(!sample.delta){
                    CameraConnection cc = scene.camera->connect_light_path(sampler, intersection);


                    //TODO: check if copy
                    LightPath new_lp = path;
                    //NOT TRANSPOSED
                    //computing p_Pt at second to last vertex, because now we
                    // can compute the pdf
                    {

                    new_lp.at(new_lp.size()-2).p_pt_acc = compute_area_pdf(new_lp.at(new_lp.size()-2), new_lp.at(new_lp.size()-1), cc.pos,false);

                    }
                    //TODO: check
                    new_lp.at(new_lp.size()-1).p_pt_acc = 1.0/cc.factor; //todo: add factor

                    Ray shadow_ray(intersection.hitpoint,
                        normalized(cc.pos-intersection.hitpoint),
                            length(cc.pos-intersection.hitpoint)-2.0*EPS);

                    shadow_ray.o = shadow_ray.o+shadow_ray.d*EPS;

                    if(bsdf->non_zero(intersection,ray.d*(-1.0),shadow_ray.d) && cc.i >-1){
                    if(!scene.primitive->intersect_any(shadow_ray)){


                        Vec3 color =  mul*bsdf->eval(
                            normalized(ray.d*(-1.0)),
                            normalized(shadow_ray.d),
                            intersection
                        ) * correcting_factor * fabs(dot(intersection.g_normal,shadow_ray.d)) / fabs(dot(intersection.normal,shadow_ray.d));
                        color = color*factor;
                        color = color*cc.factor;

                        //note: cameras have lens with zero area

                        float pdf_sum = 0.0;
                        LightPath nlp  = pdfSum(LightPath(), new_lp,pdf_sum);
                        float w_mis;
                        if(nlp.at(0).p_lt_acc>1e10){
                            w_mis = 1.0f;
                        } else {
                            w_mis = nlp.at(0).p_lt_acc/pdf_sum;
                        }

                        if(w_mis<0.0){
                            std::cout<<"NOOOO"<<std::endl;
                        }
                        if(w_mis>1.0){
                            std::cout<<"NOO222"<<std::endl;
                        }
                        //MIS
                        color = color*w_mis;
                        if(std::isnan(color.x) ||std::isnan(color.y)  || std::isnan(color.z)){
                            std::cout<<"What :("<<std::endl;
                        }


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


                {

                float r2 = lensqr(path.back().it.hitpoint - intersection.hitpoint);
                float c = fabs(dot(path.back().it.g_normal,ray.d));

                float p_pdf =  bsdf->pdf(sample.wi, ray.d*(-1.0),intersection);
                Vec3 eval_pt = bsdf->eval(sample.wi, ray.d*(-1.0),intersection);
                p_pdf*=russian_roulette(eval_pt/p_pdf);

                //TODO: no *=?
                path.back().p_pt_acc =p_pdf*r2/c;

                }

                //todo: russian roulletes only for 2 boucne or something
                if(i>-1){
                    float rr = russian_roulette(eval/p);
                    if(sampler.sample() > rr) {
                        break;
                    }
                    p = p*rr;
                }
                old_p = p;

                mul = mul*eval/p* correcting_factor * fabs(dot(intersection.g_normal,sample.wi)) / fabs(dot(intersection.normal,sample.wi));

                ray = Ray(intersection.hitpoint, sample_direction);
                ray.o = ray.o+ray.d*EPS;
                prev_intersection = intersection;
            }
        }

        return path;


    }

    float get_p_from_camera(Intersection& it){
        //todo
    }

    Vec3 integrate(const Scene& scene, const Vec2& coords, Sampler& sampler, int n_sample=0) const {
        float angular_pdf = 0.0;

        auto light_path = lightpaths.at(n_sample);


        Ray ray = scene.camera->get_ray(coords,angular_pdf);
        Intersection intersection;
        Vec3 col(0.0);
        Vec3 mul(1.0);
        bool sampled_delta = true;
        //incorrect, unless pixels are 1:1
        float p = angular_pdf;
        LightPath path;

        LightPathVertex lv;
        lv.it.hitpoint = ray.o;

        lv.p_pt_acc = 1.0;// single surface point possible: origin
        lv.p_lt_acc = 1.0; //every light path will deterministically connect to center

        //TODO: init camera vertex
        path.push_back(lv);

        //TODO: convert the probability factor into a

        for(int i = 0; i<MAX_PATH_LENGTH; i++){
            bool intersected = scene.primitive->intersect(ray,&intersection);



            if(!intersected){
                //col = col + mul*render_sky(ray)*compute_pt_direct_factor();
                break;
            } else {
                auto bsdf = intersection.get_bxdf();
                LightPathVertex lv;

                float r2= lensqr(intersection.hitpoint - path.back().it.hitpoint);

                Vec3 dir = ray.d;
                float cosTheta = fabs(dot(dir,intersection.g_normal));
                lv.p_pt_acc = compute_area_pdf(p,r2,cosTheta);

                if(lv.p_pt_acc>10000000.0){
                    std::cout<<"whaaaaa"<<std::endl;
                }

                lv.it = intersection;
                lv.bxdf = bsdf.get();
                path.push_back(lv);

                if(bsdf->is_emitter()){
                    if(sampled_delta){
                        //TODO: handle in bdpt

                        //scene.light = dynamic_cast
                        //float p_bsdf_nee = scene.light->pdf(static_cast<std::shared_ptr<EmissiveBxDF>>(bsdf))
                        //float p_bsdf_nee  = scene.light->pdf()*(r*r)/bsdf_cos;
                        Vec3 rad = bsdf->emit(ray.d*(-1.0),intersection);
                        col = col + mul*rad;
                    } else {
                        auto lid = intersection.primitive->get_light()->get_group();
                        float p_nee = scene.light->pdf(lid.first,lid.second, intersection.texture_coords,ray.o);

                        Vec3 rad = bsdf->emit(ray.d*(-1.0),intersection);
                        path.back().p_lt_acc = p_nee;
                        float sum = 0.0;
                        LightPath lp = pdfSum(path, LightPath(), sum);
                        float pdf = lp.back().p_pt_acc;
                        float w_mis = pdf/sum;

                        {
                            float r = length(intersection.hitpoint-path.at(path.size()-2).it.hitpoint);

                            //auto lid = lv.light->get_group();
                            float cos_theta = fabs(dot(intersection.g_normal,ray.d));
                            float p_nee = scene.light->pdf(lid.first,lid.second, intersection.texture_coords,ray.o);
                            float p_pt = p;
                            p_pt *=cos_theta/(r*r);
                            float w_mis2 = p_pt/(p_nee+ p_pt);


                        }
                        if(w_mis<0.0){
                            std::cout<<"NOOOO"<<std::endl;
                        }
                        if(w_mis>1.0){
                            std::cout<<"NOO222"<<std::endl;
                        }
                        col = col + mul*rad*w_mis;
                        if(std::isnan(col.x)||std::isnan(col.y)|| std::isnan(col.z)){
                            std::cout<<"whatttt2"<<std::endl;
                        }
                    }
                    break; //emitters do not reflect!
                }
                sampled_delta = false;



                DirectionalSample sample = bsdf->sample(sampler, ray.d*(-1.0), intersection);

                Vec3 sample_direction = sample.wi;
                p = sample.pdf;


                /* BDPTConnect */

                if(!sample.delta){

                    LightPath n_lp;

                    for(int lv_idx = 0; lv_idx<1; lv_idx++){
                        auto lv = light_path.at(lv_idx);
                        n_lp.push_back(lv);
                        LightSample ls;
                        ls.delta = false;
                        ls.infinite = false;
                        ls.intensity = lv.radiance;
                        ls.pdf = 1.0;
                        ls.position = lv.it.hitpoint;

                        Vec3 s_dir = ls.position-intersection.hitpoint;
                        float r = length(s_dir);
                        s_dir = s_dir/r;
                        Ray shadow_ray(intersection.hitpoint+s_dir*EPS,s_dir,r-2.0*EPS);

                        if(bsdf->non_zero(intersection,ray.d*(-1.0),shadow_ray.d)){
                        if(!scene.primitive->intersect_any(shadow_ray)){
                            Vec3 eval1 = bsdf->eval(
                                ray.d*(-1.0),
                                shadow_ray.d,
                                intersection
                            );
                            Vec3 rad(0.0);
                            if(lv_idx>0){
                                Vec3 eval2 = bsdf->eval(
                                    shadow_ray.d*(-1.0),
                                    normalized(light_path.at(lv_idx-1).it.hitpoint- light_path.at(lv_idx).it.hitpoint),
                                    lv.it
                                );
                                rad = lv.radiance*eval2;

                            } else {
                                rad = lv.light->get_emission(s_dir*(-1.0), &lv.it);//*fabs(dot(shadow_ray.d,lv.it.g_normal));
                            }

                            //todo: divide by pdf
                            float sum=0.0;
                            auto lp = pdfSum(path,n_lp,sum);
                            float pdf = lp.at(i+1).p_pt_acc*lp.at(i+2).p_lt_acc;
                            float w_mis = pdf/sum;
                            if(pdf>1e10){
                                w_mis = 1.0;
                            }

                            {
                                auto lid = lv.light->get_group();
                                float cos_theta = fabs(dot(lv.it.g_normal,shadow_ray.d));
                                float p_nee = scene.light->pdf(lid.first,lid.second, lv.it.texture_coords,ray.o);
                                float p_pt = intersection.get_bxdf()->pdf(ray.d*(-1.0f),shadow_ray.d, intersection);
                                float rr = russian_roulette(eval1/p_pt);
                                p_pt *=cos_theta/(r*r);
                                p_pt *= rr;
                                float w_mis2 = p_nee/(p_nee+ p_pt);
                            }

                            //convert p_lt_acc to hemisphere p
                            //eval already includes the cosine term;
                            float p_nee = lp.at(i+2).p_lt_acc*r*r;

                            if(w_mis<0.0){
                            std::cout<<"NOOOO"<<std::endl;
                            }
                            if(w_mis>1.0){
                                std::cout<<"NOO222"<<std::endl;
                            }
                            col =col+ mul*rad*eval1*w_mis/p_nee;
                            if(std::isnan(col.x)||std::isnan(col.y)|| std::isnan(col.z)){
                            std::cout<<"whatttt3"<<std::endl;
                            }




                        }
                        }


                    }


                } else {
                    sampled_delta = true;
                }

                /* MIS */

                if(!bsdf->non_zero(intersection,ray.d*(-1.0),sample_direction)){
                    break;
                }

                Vec3 eval = bsdf->eval(ray.d*(-1.0), sample_direction, intersection);

                //Todo: nonlinear rr
                if(i>-1){
                    float rr = russian_roulette(eval/p);
                    if(sampler.sample() > rr) {
                        break;
                    }
                    p = p*rr;
                }
                mul = mul*eval/p;
                {
                    //TODO111
                }



                ray = Ray(intersection.hitpoint+sample_direction*EPS, sample_direction);
            }
        }
        return col;
    }

public:
    VCM(const RenderSettings& rs): samples(rs.spp) {}

    static const char* name(){
        return "VCM";
    }

    void render(Scene& sc, std::string filename) {
        RgbImage* image= &(sc.camera->get_image());
        Vec2i res = image->res;
        std::vector<std::vector<LightPath>> lps(omp_get_max_threads());
        for(int k =0; k<samples*samples; k++){
            std::cout<<"rendering sample "<<k<<std::endl;

            #pragma omp parallel for
            for(int i = 0; i<lps.size();i++){
                lps.at(i) = std::vector<LightPath>();
            }
            #pragma omp parallel for
            for(int i = 0; i<res.x; i++){
                RandomSampler s(k*res.x+i);
                for(int j = 0; j<res.y; j++){
                    Vec3 acc(0.0);
                    auto lp = emit_light_path(sc,s,1.0/(samples*samples));
                    lps.at(omp_get_thread_num()).push_back(lp);

                }
            }


            //store all light paths
            lightpaths = std::vector<LightPath>();
            for(int i =0; i<lps.size(); i++){
                for(int j =0; j<lps.at(i).size();j++){
                    lightpaths.push_back(lps.at(i).at(j));
                }
            }


            //kdtree.build();

          //  #pragma omp parallel for
            for(int i = 0; i<res.x; i++){


                RandomSampler s(i+k*res.x);
                for(int j = 0; j<res.y; j++){
                    Vec3 acc(0.0);
                    float r1 = s.sample();
                    float r2 = s.sample();

                    Vector2<float> coords( ((float(i)+r1)/float(res.x))*2.0-1.0, -(((float(j)+r2)/float(res.y))*2.0-1.0));
                    acc = integrate(sc, coords, s,(i*res.y+j));
                    acc = acc/(float(samples*samples));
                    acc=acc+image->get_pixel_v(i,j);
                    image->put_pixel(i,j,acc);

                }
            }
        }


        image->save(filename);
    }
};
}
#endif