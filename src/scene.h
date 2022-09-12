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
#include<sstream>
namespace cpppt{

struct SceneData{
    std::vector<std::shared_ptr<Mesh>> meshes;
    std::vector<std::shared_ptr<Material>> materials;
};

class Scene {

public:
    std::shared_ptr<Camera> camera;
    std::shared_ptr<BVH> primitive;
    std::shared_ptr<LightGroup> light;
    /*

    static int process_material(Scene* s, std::string directory, std::string filename){

        std::ifstream ifs(directory+filename);
        std::string albedo_map;
        std::string normal_map;
        std::string specular_map;
        std::string roughness_map;
        std::string metalness_map;
        std::string transparency_map; //transparency
        //todo: every material is disney! for now, only albedo
        ifs>>albedo_map;
        std::shared_ptr<Texture> tex(get_texture(directory,albedo_map));


    }

    static Scene from_file(std::string directory, std::string filename){
        std::string line;
        std::ifstream ifs(directory+filename);
        Scene sc;
        sc.primitive = std::make_shared<BVH>();
        sc.light = std::make_shared<LightGroup>();

        while(std::getline(ifs,line)){
            std::string word;
            ifs>>word;
            if(word=="obj"){
                std::string filename;
                ifs>>filename;
                int id_mat = process_material(&sc,directory,filename);
                ifs>>filename;
                int mesh = process_mesh(&sc,directory,filename,id_mat);
            }else if(word=="light"){
                //only pointlights for now
                ifs>>word;
                if(word =="point"){

                    Vec3 pos;
                    Vec3 color;
                    float intensity;
                    ifs>>pos.x>>pos.y>>pos.z>>color.x>>color.y>>color.z>>intensity;
                    auto a = std::make_shared<PointLight>(pos,color,intensity);
                    sc.light->add(a);

                }
            } else if (word=="camera"){
                Vec2i res;
                float tan_fovy;
                Vec3 origin;
                Vec3 forward;
                Vec3 up;
                ifs>>res.x>>res.y>>tan_fovy>>origin.x>>origin.y>>origin.z>>forward.x>>forward.y>>forward.z>>up.x>>up.y>>up.z;

                auto c = std::make_shared<CameraPerspective>(res,tan_fovy,origin,forward,up);
            }
        }
    }*/
};

}
#endif