#ifndef CPPPT_RENDERER
#define CPPPT_RENDERER


#include <camera/camera.h>
#include <string>
#include <scene.h>
#include <functional>
#include <optional>
#include <omp.h>


namespace cpppt{


struct RenderSettings{
    enum RT_ENGINE{
        BVH,
        EMBREE
    };
    int renderer = 0;
    int spp = 1;
    int resX = 256;
    int resY = 256;
    std::string scene_name = "";
    std::string output_name = "";
    int max_path_length = 32;
    RT_ENGINE rt_engine=RT_ENGINE::EMBREE;
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




    //this should be everywhere but here
    static int toHemiDirection(Vec2 coords,const Intersection& is, Ray& r){
        float l = length(coords);
        if(l>=1.0){
            return -1;
        }

        float z = sqrt(1.0-l*l);
        r.o = is.hitpoint+is.normal*EPS;
        r.d = {coords.x, coords.y, z};

        //orthogonal(is.normal,&x,&y,&z);
        //Mat3 coords(y,z,x);
        Mat3 coords3(is.tangent,is.bitangent,is.normal);

        //TODO SHOULD I TRANSPOSE COORDS
        r.d = coords3*r.d;

        return 0;

    }

    //visualises some function around an hemisphere
    //this should be: a specific renderer and a camera
    template
    <typename Fn>
    static void getHemiRender(
        Fn fn, std::shared_ptr<RgbImage> image,
        int spp_sqrt=1,
        std::optional<std::function<void()>>&& begin = std::nullopt,
        std::optional<std::function<void()>>&& end = std::nullopt
    )
    {

        if(begin.has_value()){
            begin.value()();
        }

        #pragma omp parallel for
        for(int i = 0; i<image->res.x; i++){
            RandomSampler s(i);
            for(int j = 0; j<image->res.y; j++){
                Vec3 acc(0.0);
                for(int k = 0; k<spp_sqrt; k++){
                    for(int l=0; l<spp_sqrt; l++){
                        Vec2 disp = {k+s.sample(), l+s.sample()};
                        disp = disp*(1.0/spp_sqrt);
                        Vec2 coords = {
                                (float(i)+disp.x)/float(image->res.x) * 2.0f -1.0f,
                                (float(j)+disp.y)/float(image->res.y) * 2.0f -1.0f
                        };

                        acc = acc + fn(coords,s)/(spp_sqrt*spp_sqrt);

                    }

                }

                *(image->get_pixel(i,j)) = acc;


            }
        }

        if(end.has_value()){
            end.value()();
        }
    }

/*
    template
    <typename Fn>
    static void getHemiSplat(
        Fn fn, RgbImage& image,
        int spp_sqrt=1,
        std::optional<std::function<void()>>&& begin = std::nullopt,
        std::optional<std::function<void()>>&& end = std::nullopt
    )
    {

        if(begin.has_value()){
            begin.value()();
        }

        #pragma omp parallel for
        for(int i = 0; i<image.res.x; i++){
            RandomSampler s(i);
            for(int j = 0; j<image.res.y; j++){
                for(int k = 0; k<spp_sqrt; k++){
                    for(int l=0; l<spp_sqrt; l++){
                        coords = fn(col,s);
                        col = col/(spp_sqrt*spp_sqrt);
                        Vec3* pix = *(image.get_pixel(i,j));

                        #pragma omp atomic
                        pix.x += col.x;
                        #pragma omp atomic
                        pix.y += col.y;
                        #pragma omp atomic
                        pix.z += col.z;

                    }

                }


            }
        }

        if(end.has_value()){
            end.value()();
        }
    }
    */




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