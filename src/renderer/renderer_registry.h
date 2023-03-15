#ifndef CPPPT_RENDERER_REGISTRY
#define CPPPT_RENDERER_REGISTRY


#include <camera/camera.h>
#include <string>
#include <scene.h>
#include <unordered_map>
#include <renderer/renderer.h>
#include <renderer/pathtracer.h>
#include <renderer/pathtracerMIS.h>
#include <renderer/lighttracer.h>
#include <renderer/debug_renderer.h>
#include <renderer/pathtracerMLT.h>
#include <renderer/lighttracerMLT.h>
#include <renderer/photonmapping.h>
#include <renderer/ppg.h>
#include <renderer/VCM.h>



namespace cpppt{


class AbstractFactory{
public:
    virtual std::shared_ptr<Renderer> create(const RenderSettings& rs) const = 0;
};

template
<typename T>
class RendererFactory: public AbstractFactory{
public:
    std::shared_ptr<Renderer> create(const RenderSettings& rs) const override {
        return std::make_shared<T>(rs);
    }
};


struct RendererRegistry{
    std::vector<std::unique_ptr<AbstractFactory>> factories;
    std::vector<const char*> names;
    std::unordered_map<const char*, size_t> id;


    template
    <typename T>
    void reg(){
        names.push_back(T::name());
        factories.push_back(std::make_unique<RendererFactory<T>>());
        id.insert({T::name(), names.size()-1});
    }

    std::shared_ptr<Renderer> create(const char* name, const RenderSettings& rs){
        return create(id[name], rs);
    }
    std::shared_ptr<Renderer> create(int i,const RenderSettings& rs){
        return factories.at(i)->create(rs);
    }

    static std::shared_ptr<RendererRegistry> get(){

        return RendererRegistry::instance;

    }

    RendererRegistry(){
        reg<Pathtracer>();
        reg<PathtracerMIS>();
        reg<Lighttracer>();
        reg<DebugRenderer>();
        reg<PathtracerMLT>();
        reg<LighttracerMLT>();
        reg<PhotonMapping>();
        reg<PPG>();
    }

    private:

        static inline std::shared_ptr<RendererRegistry> instance
            =std::make_shared<RendererRegistry>();

};
}
#endif