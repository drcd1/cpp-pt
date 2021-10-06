#ifndef CPPPT_SCENE
#define CPPPT_SCENE

#include <camera/camera.h>
#include <primitive/primitive.h>

namespace cpppt{
struct Scene {
    Camera* camera;
    Primitive* primitive;
};
}
#endif