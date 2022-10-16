#ifndef CPPPT_RENDERER
#define CPPPT_RENDERER


#include <camera/camera.h>
#include <string>
#include <scene.h>

namespace cpppt{

struct RenderSettings{
    enum RendererType {PATHTRACER, PATHTRACER_MIS, LIGHTTRACER, DEBUG, PATHTRACER_MLT, LIGHTTRACER_MLT};

    RendererType renderer = RendererType::DEBUG;
    int spp = 1;
    int resX = 256;
    int resY = 256;
    std::string scene_name = "";
    std::string output_name = "";
};

class Renderer{
public:
    virtual void render(Scene& sc, std::string filename) const = 0;
};
}

#endif