#ifndef CPPPT_RENDERER
#define CPPPT_RENDERER


#include <camera/camera.h>
#include <string>
#include <scene.h>
namespace cpppt{
class Renderer{
    virtual void render(Scene& sc, std::string filename) const = 0;
};
}

#endif