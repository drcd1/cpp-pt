#ifndef CPPPT_RENDERER
#define CPPPT_RENDERER


#include <camera/camera.h>
#include <string>
#include <scene.h>
#include <functional>


namespace cpppt{


struct RenderSettings{
    int renderer = 0;
    int spp = 1;
    int resX = 256;
    int resY = 256;
    std::string scene_name = "";
    std::string output_name = "";
};

struct RenderProgress{
    float percentage;
    int current_task;
    int total_tasks;
    char task_description[128]="\0";
};

class RenderMessenger{
public:
    virtual void eval() const = 0;
};


class Renderer{
private:
    RenderMessenger* messenger=nullptr;

public:
    virtual void render(Scene& sc, std::string filename) = 0;
    void renderW(Scene& sc, std::string filename){
        render(sc,filename);
        finish();
    };

    virtual void getRayRadiance(const Ray& r, const Scene& sc){

    }

    

    //this should be everywhere but here
    static int toHemiDirection(Vec2 coords,const Intersection& is, Ray& r){
        float l = coords.length();
        if(coords.length>=1.0){
            return -1;
        }
        
        float z = sqrt(1.0-coords.lensqr());
        r.o = is.hitpoint;
        r.d = {coords.x, coords.y, z};

        //orthogonal(is.normal,&x,&y,&z);
        //Mat3 coords(y,z,x);
        Mat3 coords(is.tangent,is.bitangent,is.normal);

        //TODO SHOULD I TRANSPOSE COORDS
        r.d = coords*r.d; 

        return 0;

    }

    //visualises some function around an hemisphere
    //this should be: a specific renderer and a camera
    static getHemiRender(
        std::function<Vec3(const Vec2& coords)>& fn, RgbImage& image, 
        int spp_sqrt, 
        std::optional<std::function<void()>>&& begin = std::nullopt,
        std::optional<std::function<void()>>&& end = std::nullopt   
    )
    {

        if(begin.has_value()){
            begin.value()();
        }

        #pragma omp_parallel_for
        for(int i = 0; i<image.width(); i++){
            Sampler s(i);
            for(int j = 0; j<image.height(); j++){
                Vec3 acc;
                for(int k = 0; k<spp_sqrt; k++)
                    for(int l=0; l<spp_sqrt; l++){
                        Vec2 disp = {k+s.sample(), j+s.sample()}/spp_sqrt;
                        Vec2 coords = {
                                (i+disp.x)/image.width * 2.0 -1.0,
                                (j+disp.y)/image.width * 2.0 -1.0
                        }

                        acc +=fn(coords)/(spp_sqrt*spp_sqrt);
                        
                    }

                }

                image.get_pixel(i,j) = acc;   


            }
        } 
        
        if(end.has_value()){
            end.value()();
        }       
    }




    virtual void updateProgress(RenderProgress& progress){
        progress.percentage = 0.0;
        progress.current_task = 0;
        progress.total_tasks = 0;
        sprintf(progress.task_description,"Progress Unavailable For This Renderer");
    }

    void setRendererMessenger(RenderMessenger* m){
        messenger = m;
    }

    void finish(){
        if(messenger)
           messenger->eval();
    }
};
}

#endif