#include <camera/camera_perspective.h>
#include <renderer/dummy_renderer.h>
#include <renderer/debug_renderer.h>
#include <renderer/Pathtracer.h>
#include <renderer/Lighttracer.h>
#include <renderer/pathtracerMLTstrat1.h>
#include <primitive/simple_group.h>
#include <primitive/bvh.h>
#include <shape/mesh.h>
#include <shape/triangle.h>
#include <primitive/primitive_leaf.h>
#include <bxdf/diffuse_bxdf.h>
#include <bxdf/emissive_bxdf.h>
#include <texture/constant_texture.h>
#include <bxdf/mirror_bxdf.h>
#include <bxdf/refraction_bxdf.h>

#include <light/light.h>
#include <light/point_light.h>
#include <light/shape_light.h>
#include <light/light_group.h>
#include <loader/loader.h>
#include <shape/shape.h>
#include <scene.h>
#include <chrono>


using namespace cpppt;

int main(int argc, char** argv){
    try{
    if(argc<4){
        std::cout<<"Insuficient Arguments"<<std::endl;
        return 0;
    }

    RenderSettings rs;
    Scene s;
    SceneData sd;
    Loader::load_render_settings(&rs,argv[1]);
    Loader::load_scene(&s,&sd,&rs,argv[2]);

    std::unique_ptr<Renderer> renderer;

    if(rs.renderer == RenderSettings::RendererType::PATHTRACER){
        renderer = std::make_unique<Pathtracer>(rs);
    } else if(rs.renderer == RenderSettings::RendererType::LIGHTTRACER){
        renderer = std::make_unique<Lighttracer>(rs);
    } else if(rs.renderer == RenderSettings::RendererType::DEBUG){
        renderer = std::make_unique<DebugRenderer>();
    }  else if(rs.renderer == RenderSettings::RendererType::PATHTRACER_MLT){
        renderer = std::make_unique<PathtracerMLT>(rs);
    }

    //DummyRenderer renderer;


    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    auto t1 = high_resolution_clock::now();

    renderer->render(s, argv[3]);

    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;

    std::cout << ms_double.count()*0.001 << "s\n";
    } catch(std::runtime_error e){
        std::cout<<e.what()<<std::endl;
    }
}
