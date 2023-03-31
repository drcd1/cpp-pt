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

        lp.push_back(v1);
        lp.push_back(v2);
        lp.push_back(v3);

        }




    }


    float compute_area_pdf(const LightPathVertex& v, const LightPathVertex& v_from, const Vec3& from_from, bool transposed){
        //TODO: account for transposed
        Vec3 wo = normalized(from_from-v_from.pt);
        Vec3 wi = v.pt-v_from.pt;

        float r= length(wi);
        wi = wi/r;

        float cosTheta = fabs(dot(wi,v.it.g_normal));

        float p = v_from.bxdf.pdf(wo,wi,v_from.it);
        float rr = russian_roulette(v_from.bxdf.eval(wo,wi,v_from.it));

        return compte_area_pdf(p,r*r,cosTheta);
    }
    float compute_area_pdf(const LightPathVertex& v, const LightPathVertex& v_from, float solid_angle_pdf){
        //TODO: account for transposed
        Vec3 wi = v.pt-v_from.pt;

        float r= length(wi);
        wi = wi/r;
        float cosTheta = fabs(dot(wi,v.it.g_normal));

        return compte_area_pdf(solid_angle_pdf,r*r,cosTheta);
    }

    inline float compute_area_pdf_r(float solid_angle_pdf, float r, float cosTheta){
        return compute_area_pdf(solid_angle_pdf,r*r,cosTheta);
    }
    
    inline float compute_area_pdf(float solid_angle_pdf, float rSQ, float cosTheta){
        return solid_angle_pdf*cosTheta/rSQ;
    }


    //assumes every vertex has the conditional
    //probabilities expressed in terms of area
    void accPdf(LightPath& lp){
        for(int i = 1; i<lp.size(); i++){
            lp.at(i).p_pt_acc *= lp.at(i-1)*p_pt_acc;
            lp.at(lp.size()-i-1).p_lt_acc *= lp.at(lp.size()-i);
        }
    }


    LightPath pdfSum(const LightPath& from_camera, const LightPath& from_light, float& pdfSum){
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
        float wo = normalized(ll.previous->pt - ll->pt);
        float wi = dir;
        
        float p = ll.bxdf->pdf(wo,wi,ll.it);
        p*= russian_roullette(eval(ll.bxdf->eval(wo,wi,ll.it)));
        //lc.p_lt_acc = ll.p_lt_acc*fabs(dot(dir,lc.it.g_normal))/(r*r);
        lc.p_lt_acc =  p*fabs(dot(dir,lc.it.g_normal))/(r*r);

        }
        //NOT TRANSPOSED
        {        
        float wi = dir;
        float wo = normalized(lc.previous->pt - lc->pt);

        float p = lc.bxdf->pdf(wi,wo,lc.it);
        p*= russian_roullette(eval(ll.bxdf->eval(wo,wi,ll.it)));
        
        ll.pt_acc = p*fabs(dot(dir,ll.it.g_normal))/(r*r);
        }

        //Now all the paths probabilities are computed, and we can just sum


        //todo: okay, so this is way to0 many copies
        LightPath l = from_camera;
        l.insert(l.end(),from_light.rbegin(),from_light.rend());
        accPdf(l);

        pdfSum = computePdfSumFullPath(l);
        return l;
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
        lpv.bxdf = nullptr; //todo
        //lpv.intersection = nullptr;//todo

        path.push_back(lpv);


        Intersection intersection;
        Intersection prev_intersection;
        Vec3 mul = lps.radiance/lps.pdf;

        RgbImage& image = scene.camera->get_image();
        float old_p = lps.angle_pdf;
        for(int i = 0; i<MAX_PATH_LENGTH; i++){
            bool intersected = scene.primitive->intersect(ray,&intersection);
            if(!intersected){
                return path;
            } else {
                auto bsdf = intersection.get_bxdf();
                if(bsdf->is_emitter()){
                    break;
                }

                LightPathVertex lpv;
                lpv.it = intersection;

                
                {
                
                lpv.p_lt_acc = compute_area_pdf(old_p,lpv,path.back());
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


                    //TODO: check if copy
                    LightPath new_lp = path;
                    //NOT TRANSPOSED
                    //computing p_Pt at second to last vertex, because now we
                    // can compute the pdf
                    {
                    
                    lp.at(lp.size()-2).p_pt_acc = compute_area_pdf(lp.at(lp.size()-2), lp.at(lp.size()-1), cc.pos,false);

                    }   
                    //TODO: check
                    lp.at(lp.size()-1).p_pt_acc = 1.0/cc.factor;

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
                      /*  for(int k = path.size()-1; k>0; k--){
                            //when computing the weights, k is the last vertex generated by pathtracing
                            //k = 0 is pathtracing, not counted in this loop
                            //k = 1 is nee
                            //k = path.size() is lighttracing, not counted in this loop
                            p_forward*=path.at(k).p_pt_acc;
                            denom += p_forward*path.at(k).p_pt_acc*connection_prob*path.at(k).p_lt_acc;
                        }
                        denom += p_forward*path.at(0).p_pt_acc;
                        float w = p/denom;
                        
                        */

                        LightPath lp = pd




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


                {

                float r2 = lensqr(path.back().pt - intersection.hitpoint);
                float c = fabs(dot(path.back().g_normal,ray.d));
                //TODO: russian roulette prob
                
                path.back().p_pt_acc*= bsdf->pdf(sample.wi, ray.d*(-1.0),intersection)*r2/c;

                }

                //todo: russian roulletes only for 2 boucne or something
                if(i>-1){
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

    float get_p_from_camera(Intersection& it){
        //todo
    }

    Vec3 integrate(const Scene& scene, const Vec2& coords, Sampler& sampler) const {
        Ray ray = scene.camera->get_ray(coords);
        Intersection intersection;
        Vec3 col(0.0);
        Vec3 mul(1.0);
        bool sampled_delta = true;
        //incorrect, unless pixels are 1:1
        float p = 1.0;
        LightPath path;


        std::vector<LightPathVertex> lc;
        LightPathVertex lv;

        lv.p_pt_acc = 1.0;// single surface point possible: origin
        lv.p_lt_acc = 1.0; //every light path will deterministically connect to center

        //TODO: init camera vertex
        lc.push_back(lv);

        //TODO: convert the probability factor into a 

        for(int i = 0; i<MAX_PATH_LENGTH; i++){
            bool intersected = scene.primitive->intersect(ray,&intersection);



            if(!intersected){
                //col = col + mul*render_sky(ray)*compute_pt_direct_factor();
                break;
            } else {
                auto bsdf = intersection.get_bxdf();
                LightPathVertex lv;

                Vec3 r2= (intersection.hitpoint - lc.back().it.hitpoint).lensqr();
                if(i>0){
                    lv.p_pt_acc = p*fabs(dot(dir,intersection.g_normal))/r2;
                } else {
                    lv.p_pt_acc = get_p_from_camera(intersection);
                }

                lv.it = intersection;
                lv.bxdf = bsdf;
                lc.push_back(lv)

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
                        lc.back().p_lt_acc = p_nee;
                        float sum = 0.0;
                        LightPath lp = pdfSum(lc, LightPath(), &sum);
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
                    for(lv: path)


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
        std::vector<std::vector<LightPath>> lps(omp_get_max_threads);
        for(int k =0; k<samples*samples; k++){
        #pragma omp parallel for
        for(int i = 0; i<res.x; i++){



            RandomSampler s(k*res.x+i);
            for(int j = 0; j<res.y; j++){
                Vec3 acc(0.0);
                auto lp = emit_light_path(sc,s,1.0/(res.x*res.y));
                lps.at(omp.get_thread_num()).push_back(lp);

            }
        }

        //store all light paths

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