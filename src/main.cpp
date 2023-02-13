#include <camera/camera_perspective.h>
#include <renderer/dummy_renderer.h>
#include <renderer/debug_renderer.h>
#include <renderer/Pathtracer.h>
#include <renderer/Lighttracer.h>
#include <renderer/pathtracerMLT.h>
#include <renderer/pathtracerMIS.h>
#include <renderer/lighttracerMLT.h>
#include <renderer/photonmapping.h>
#include <renderer/pure_pt.h>
#include <renderer/ppg.h>
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
#include <time.h>
#include <algorithms/sdtree.h>

#include <light/light.h>
#include <light/point_light.h>
#include <light/shape_light.h>
#include <light/light_group.h>
#include <loader/loader.h>
#include <shape/shape.h>
#include <scene.h>
#include <chrono>

//NOTE: sometimes, sample is not zero, but prob is zero?!
//I think we are dealing with this in the renderer by discarding invalid samples
//but i need to check if this only happens for invalid samples
// if so, why nans?
using namespace cpppt;
/*
void runDTreeTest(){

    Bounds<Vec3> bounds;
    bounds.min = Vec3(0.0);
    bounds.max = Vec3(1.0);
    SDTree sdt(bounds);

    sdt.splat(0.0)


}
*/
void runValidation(){
    Intersection it;
    it.normal = Vec3(0,0,1);
    it.tangent = Vec3(1,0,0);
    it.bitangent = Vec3(0,1,0);
    it.g_normal= it.normal;
    it.hitpoint = Vec3(0,0,0);
    it.computed_normal = true;
    it.texture_coords = Vec3(0.0);
    it.material = std::make_shared<StandardMaterial>(
      std::make_shared<ConstantTexture<Vec3>>(Vec3(0.8f,0.9f,0.8f)),
      std::make_shared<ConstantTexture<Vec3>>(Vec3(0.5f,0.5f,1.0f)),
      std::make_shared<ConstantTexture<float>>(0.5f),
      std::make_shared<ConstantTexture<float>>(0.5f),
      std::make_shared<ConstantTexture<float>>(0.0f),
      std::make_shared<ConstantTexture<float>>(1.23f),
      std::make_shared<ConstantTexture<float>>(0.0f),
      std::make_shared<ConstantTexture<float>>(0.0f)
    );

    auto bxdf = it.get_bxdf();
    Vec3 eye = normalized(Vec3(-1.0,1.0,0.3));

    RandomSampler s(time(NULL));
    DirectionalSample sample = bxdf->sample(s,eye,it);
    Vec3 light = sample.wi;
    float p = sample.pdf;

    std::cout<<"Light sampled: "<<light.x<<" "<<light.y<<" "<<light.z<<std::endl;
    std::cout<<"Prob: "<<p<<std::endl;
    std::cout<<"pdf: "<<bxdf->pdf(eye,light,it)<<std::endl;
    Vec3 eval = bxdf->eval(eye,light,it);
    std::cout<<"eval: "<<eval.x<<" "<<eval.y<<" "<<eval.z<<std::endl;


}

int main(int argc, char** argv){
    try{
    if(argc<2){
        std::cout<<"Insuficient Arguments"<<std::endl;
        return 0;
    }

    RenderSettings rs;
    Scene s;
    SceneData sd;
    Loader::load_render_settings(&rs,argv[1]);
    Loader::load_scene(&s,&sd,&rs,rs.scene_name);

    std::unique_ptr<Renderer> renderer;

    if(rs.renderer == RenderSettings::RendererType::PATHTRACER){
        renderer = std::make_unique<Pathtracer>(rs);
    } else if(rs.renderer == RenderSettings::RendererType::PATHTRACER_MIS){
        renderer = std::make_unique<PathtracerMIS>(rs);
    } else if(rs.renderer == RenderSettings::RendererType::LIGHTTRACER){
        renderer = std::make_unique<Lighttracer>(rs);
    } else if(rs.renderer == RenderSettings::RendererType::DEBUG){
        renderer = std::make_unique<DebugRenderer>();
    }  else if(rs.renderer == RenderSettings::RendererType::PATHTRACER_MLT){
        renderer = std::make_unique<PathtracerMLT>(rs);
    } else if(rs.renderer == RenderSettings::RendererType::LIGHTTRACER_MLT){
        renderer = std::make_unique<LighttracerMLT>(rs);
    } else if(rs.renderer == RenderSettings::RendererType::PHOTONMAPPING){
        renderer = std::make_unique<PhotonMapping>(rs);
    } else if(rs.renderer == RenderSettings::RendererType::PURE_PT){
        renderer = std::make_unique<PurePt>(rs);
    }else if(rs.renderer == RenderSettings::RendererType::PPG){
        renderer = std::make_unique<PPG>(rs);
    }else if(rs.renderer == RenderSettings::RendererType::VALIDATION){
        runValidation();
        return 0;
    }

    //DummyRenderer renderer;


    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    auto t1 = high_resolution_clock::now();

    renderer->render(s, rs.output_name);

    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;

    std::cout << ms_double.count()*0.001 << "s\n";
    } catch(std::runtime_error e){
        std::cout<<e.what()<<std::endl;
    }
}
