#ifndef CPPPT_RENDERER_GUI
#define CPPPT_RENDERER_GUI

#include <renderer/renderer_registry.h>
#include <thread>
#include <loader/loader.h>
#include <scene.h>
#include <sstream>
#include <imgui.h>
#include <omp.h>

namespace cpppt{

class RendererGUI;

class RenderMessengerGUI: public RenderMessenger{
    RendererGUI* rg;
public:
    RenderMessengerGUI(RendererGUI* rg):rg(rg){}

    void eval() const override ;
};



class RendererGUI{

    friend class RenderMessengerGUI;
    enum Task{
        NONE,
        LOADING,
        RENDERING
    };

    struct SceneGUI{
        char scene_filename[128] = "..\\scenes\\scene_hdr_real\\scene.ptscene";
        char output_name[128] = "scene_hdr";
        Vec2i res={256,256};
        int spp=1;

    };

    SceneGUI sceneGUI;

    int current_renderer=0;
    int old_renderer=0;
    std::shared_ptr<Renderer> renderer = nullptr;
    Scene scene;
    SceneData scene_data;
    Task current_task=Task::NONE;
    bool joined;
    std::unique_ptr<std::thread> rendering_thread= nullptr;
    RenderMessengerGUI messenger;
    RenderProgress status;
    char status_message[128];

    RenderSettings rs;
    std::shared_ptr<RendererRegistry> reg = RendererRegistry::get();

public:
    RendererGUI():messenger(this){
        omp_set_num_threads(std::max(1, omp_get_max_threads()-1));

    }
    void frame(){
        if(rendering_thread!=nullptr &&
            current_task==Task::NONE
        ){

            std::cout<<"Joining"<<std::endl;
            rendering_thread->join();
            rendering_thread = nullptr;
        }


        ImGui::Begin("CPP-PT, a pathtracer in cpp");

        ImGui::Text("This is some useful text.");

        ImGui::InputText("Scene Filename",sceneGUI.scene_filename,128);
        ImGui::InputText("Output Name",sceneGUI.output_name,128);
        if(ImGui::BeginCombo("Renderer", reg->names[current_renderer])){

            for(int i = 0; i<reg->names.size(); i++){
                const bool is_selected = (current_renderer == i);
                if(ImGui::Selectable(reg->names[i],is_selected)){
                    current_renderer = i;
                }
                // Set the initial focus when opening the combo (scrolling + keyboard navigation focus)
                if (is_selected)
                    ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }


        ImGui::DragInt("Res X", &sceneGUI.res.x);
        ImGui::DragInt("Res Y", &sceneGUI.res.y);
        ImGui::DragInt("SPP", &sceneGUI.spp);

        if(ImGui::Button("Render") && current_task == NONE){
            current_task = Task::LOADING;
            rs.resX = sceneGUI.res.x;
            rs.resY = sceneGUI.res.y;
            rs.renderer = current_renderer;
            rs.spp = sceneGUI.spp;
            rs.scene_name = std::string(sceneGUI.scene_filename);
            rs.output_name = sceneGUI.output_name;
            renderer=reg->create(current_renderer,rs);


            rendering_thread = std::make_unique<std::thread>([&](){

                std::string renderer_name = reg->names[current_renderer];

                std::stringstream ss;
                ss<<rs.output_name<<"x"<<rs.resX<<"y"
                <<rs.resY<<"spp"<<rs.spp
                <<renderer_name;

                rs.output_name = ss.str();
                //todo: only load if scene has changed
                scene = Scene();
                scene_data = SceneData();

                Loader::load_scene(&scene,&scene_data,&rs,rs.scene_name);
                current_task = Task::RENDERING;
                renderer->setRendererMessenger(&messenger);
                renderer->renderW(scene, rs.output_name);
            });

            std::cout<<"not blocked"<<std::endl;
        }

        if(renderer!=nullptr && current_task == Task::LOADING){
            ImGui::Text("Loading Scene");
        }
        if(renderer!=nullptr && current_task == Task::RENDERING){
            renderer->updateProgress(status);

            ImGui::ProgressBar(status.percentage);
            ImGui::SameLine(0.0f, ImGui::GetStyle().ItemInnerSpacing.x);
            char buf[32];
            sprintf(buf, "Task %02d/%02d",status.current_task, status.total_tasks);
            std::cout<<buf<<" "<<status.percentage<<" "<<status.task_description<<std::endl;

            ImGui::Text(buf);
            ImGui::Text(status.task_description);

        }


        ImGui::End();

    }

};

void RenderMessengerGUI::eval() const  {
    rg->current_task = RendererGUI::Task::NONE;
}


}

#endif