#ifndef CPPPT_ENVIRONMENT_LIGHT
#define CPPPT_ENVIRONMENT_LIGHT

#include <light/light.h>
#include <primitive/primitive_leaf.h>
#include <math/math.h>
#include <cmath>
#include <algorithm>
/*
For area lights
intensity is irradiance (power/area)

*/

#define SQRT_N_SAMPLES 2

namespace cpppt{
class EnvironmentLight: public Light {
private:
    int h;
    int w;
    Vec3 min;
    Vec3 max;
    std::shared_ptr<Texture> map;

    int idx(int i, int j) const {
        //todo: z-order?
        return i*h+j;
    }

    std::pair<Vec2,Vec2> coords(int idx) const {
        int i = idx/h;
        int j = idx%h;
        return {Vec2(float(i)/float(w), float(j)/float(h)),Vec2(1.0f/float(w),1.0f/float(h))};
    }

    int inv_coords(Vec2 c) const {
        c = c*Vec2(w,h);
        int i = int(c.x);
        int j = int(c.y);
        return idx(i,j);
    }

    Vec2 dir2uv(Vec3 d) const {
        float angle = atan2(d.y,d.x);
        if(angle<0.0)
            angle = angle+2.0*M_PI;
        return Vec2(angle/M_PI* 0.5,d.z*0.5+0.5);
    }

    Vec3 uv2dir(Vec2 c) const {
        c.y = c.y*2.0-1.0;
        float sinTheta = sqrt(1.0-c.y*c.y);
        return Vec3(cos(c.x*2.0*M_PI)*sinTheta,sin(c.x*2.0*M_PI)*sinTheta,c.y);
    }

    float rd(const Vec3& a) const {
        return std::max(a.x,std::max(a.y,a.z));
    }

    int chooseDirection(Sampler& s) const {
        float r = s.sample()*radiance.at(radiance.size()-1);
        int min = 0;
        int max = radiance.size();

        while(max-min>4){
            int t = (max+min)/2;
            if(radiance.at(t)>=r){
                max = t;
            } else {
                min = t;
            }
        }
        //todo: check edge cases
        for(int i = min; i<max; i++){
            if(radiance.at(i)<r && radiance.at(i+1)>=r){
                return i;
            }
        }

        throw std::runtime_error("Failed to sample hdr");
    }

    std::vector<float> radiance;



public:

    EnvironmentLight(std::shared_ptr<Texture> map,const Vec3& min, const Vec3& max,int width = 1024, int height = 512):map(map),w(width),h(height),radiance(width*height+1){

        radiance.at(0)=0.0;
        for(int i = 0; i<radiance.size()-1; i++){
            radiance.at(i+1)=0.0;
            auto c = coords(i);
            for(int k = 0; k<SQRT_N_SAMPLES; k++){
                for(int l = 0; l<SQRT_N_SAMPLES; l++){
                    Vec2 uv = c.first + c.second*Vec2((k+0.5)/SQRT_N_SAMPLES, (l+0.5)/SQRT_N_SAMPLES);
                    uv.y = uv.y*2.0-1.0;
                    float sinTheta = sqrt(1.0-uv.y*uv.y);
                    Vec3 direction =Vec3(cos(uv.x*2.0*M_PI)*sinTheta,sin(uv.x*2.0*M_PI)*sinTheta,uv.y);

                    //uv.y = (acos(uv.y*2.0-1.0)/M_PI);
                    //radiance.at(i+1) += rd(map->sample(Vec3(uv.x,uv.y,0.0)))/(SQRT_N_SAMPLES*SQRT_N_SAMPLES);
                    radiance.at(i+1)+=rd(emit(direction)/(SQRT_N_SAMPLES*SQRT_N_SAMPLES));
                }
            }
            radiance.at(i+1) += radiance.at(i);
        }




    };

    Vec3 emit(const Vec3& dir)const override  {
        return map->sample(-Vec3(atan2(dir.y,dir.x)/(2.0*M_PI)+0.5,acos(dir.z)/(M_PI),0.0));
    }

    bool is_infinite_not_delta() const override {
        return true;
    }



    //for NEE
    LightSample connect_eye_path(Sampler& s, const Intersection& from) const {
        LightSample ls;
        Intersection it;


        int c = chooseDirection(s);
        ls.pdf = (radiance.at(c+1)-radiance.at(c))/radiance.at(radiance.size()-1);
        float r1 = s.sample();
        float r2 = s.sample();


        auto p = coords(c);

        Vec2 uv = p.first + p.second*Vec2(r1,r2);


        //uv.x = (uv.x-0.5);
        uv.y=uv.y*2.0-1.0;
        float sin_theta = sqrt(1.0-uv.y*uv.y);

        Vec3 direction = Vec3(cos(uv.x*2.0*M_PI)*sin_theta,sin(uv.x*2.0*M_PI)*sin_theta,uv.y);

        ls.position = from.hitpoint + direction;

        //ls.intensity = map->sample(Vec3(uv.x,uv.y,0.0));
        //ls.pdf = 1.0/(4.0*M_PI);
        ls.intensity = emit(direction);

        ls.pdf = ls.pdf*((w*h)/(4.0*M_PI));
        ls.infinite = true;
        ls.delta = false;
        ls.normal = -direction;
        return ls;

    }

    float pdf(int l_id, const Light* parent, const Vec3& coords, const Vec3& lit) const {
        Vec3 dir = coords;
        Vec2 uv = dir2uv(dir);
        int i = inv_coords(uv);
        float pdf = (radiance.at(i+1)-radiance.at(i))/radiance.at(radiance.size()-1);
        pdf=pdf*((w*h)/(4.0*M_PI));

        return pdf;
    }

    LightPathStart sample(Sampler& s) const {
        throw std::runtime_error("Not implemented");
        LightPathStart lps;
        return lps;
    }


    bool is_delta() const {
        return false;
    }
};
}
#endif