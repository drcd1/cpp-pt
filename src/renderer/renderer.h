#ifndef CPPPT_RENDERER
#define CPPPT_RENDERER


#include <camera/camera.h>
#include <string>
#include <scene.h>



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