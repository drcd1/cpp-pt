#ifndef CPPPT_SCENE
#define CPPPT_SCENE

#include <camera/camera.h>
#include <primitive/primitive.h>
#include <primitive/bvh.h>
#include <light/light.h>
#include <light/light_group.h>
#include <shape/mesh.h>
#include <bxdf/standard_material.h>
#include <texture/texture.h>
#include <fstream>
#include <sstream>
#include <memory>
namespace cpppt{

struct SceneData{
    std::vector<std::shared_ptr<Mesh>> meshes;
    std::vector<std::shared_ptr<Material>> materials;
};

class Scene {

public:
    std::shared_ptr<Camera> camera;
    std::shared_ptr<Primitive> primitive;
    std::shared_ptr<Light> light;
};

}
#endif