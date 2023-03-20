#ifndef CPPPT_RENDERER_GUI
#define CPPPT_RENDERER_GUI

#include <renderer/renderer_registry.h>
#include <thread>
#include <memory>
#include <atomic>
#include <vickyr.h>
#include <loader/loader.h>
#include <scene.h>
#include <sstream>
#include <imgui.h>
#include <omp.h>
#include <queue>
#include <debug_tools/gui/camera.h>
#include <debug_tools/gui/clock.h>
#include <debug_tools/gui/scene.h>
#include <camera/camera_perspective.h>

#define MAX_CHARS 128

namespace cpppt{
namespace DebugTools{
class RendererGUI;

class RenderMessengerGUI: public RenderMessenger{
    RendererGUI* rg;
public:
    RenderMessengerGUI(RendererGUI* rg):rg(rg){}

    void eval() const override ;
};

struct Task{
    enum Name {
        NONE,
        LOADING,
        RENDERING
    };
    std::function<void()> task;
    Name name;
};



class RendererGUI{

    friend class RenderMessengerGUI;

    CameraGUI camera_gui;


    struct SharedVars{

        RenderSettings rs;
        std::shared_ptr<Renderer> renderer = nullptr;
        Scene scene;
        SceneData scene_data;
        bool load_render=false;

    };


    SceneGUI sceneGUI;
    SharedVars sv;

    DebugTools::Camera<DebugTools::WalkControls> camera = {
        DebugTools::WalkControls::from_vecs(
                    glm::vec3(0.0f,1.0f,0.0f) ,
                    glm::vec3(1.0f,0.0f,0.0f) ,
                    glm::vec3(0.0f,0.0f,1.0f) ,
                    glm::vec3(0.0f,-5.0f,1.0f)
                ),
                1.0
    };

    int current_renderer=0;
    int old_renderer=0;

    std::queue<Task> tasks;
    std::unique_ptr<std::thread> task_thread = nullptr;

    //helper for printer
    char status_message[MAX_CHARS];

    //real time renderer is single threaded
    VickyR::VickyR* renderer_vk;

    //registry is shared, but supposed to be const?
    //make sure it is so
    const std::shared_ptr<RendererRegistry> reg = RendererRegistry::get();

    RenderProgress status;
    //these vars are shared, but atomic
    std::atomic<Task::Name> current_task=Task::Name::NONE;
    std::atomic<bool> reloadScene = false;
    std::atomic<bool> sceneUpToDate = false;
    RenderMessengerGUI messenger;

    struct Mouse{
        double x;
        double y;

        double old_x;
        double old_y;
        bool init = false;
    };
    Clock clock;

    VickyR::AllocatedTextureImage showRenderTex;
    VkDescriptorSet renderTexId;
/*
    VickyR::TextureImage tim = renderer_vk->createTextureImage();
    renderer_vk->deleteTexture(tim);

    renderer_vk->createDummyTexture();
*/
    VkDescriptorSet tex_id = ImGui_ImplVulkan_AddTexture(
        renderer_vk->getTextureImage().sampler,
        renderer_vk->getTextureImage().view,
        renderer_vk->getTextureImage().layout

    );


    Mouse mouse;


public:
    RendererGUI(VickyR::VickyR& renderer):messenger(this),renderer_vk(&renderer){
        omp_set_num_threads(std::max(1, omp_get_max_threads()-1));
        showRenderTex = renderer_vk->getIV().imagePool.allocate();


        int w=128;
        int h=128;
        uint8_t* tmp = new uint8_t[w*h*4];
        for(int i = 0; i<h;i++){
            for(int j =0; j<w; j++){
                tmp[(i*w+j)*4] = (i%2!= j%2)?255:0;
                tmp[(i*w+j)*4+1] = 0;
                tmp[(i*w+j)*4+2] = (i%2!= j%2)?255:0;
                tmp[(i*w+j)*4+3] = 255;
            }
        }

        renderer_vk->createFullTextureImage(*showRenderTex.im, tmp,w,h);
        renderTexId = ImGui_ImplVulkan_AddTexture(
            showRenderTex.im->sampler,
            showRenderTex.im->view,
            showRenderTex.im->layout
        );

        delete[] tmp;
    }

    CameraGUI& getCameraGUI(){
        return camera_gui;
    }

    void resetMouse(){
        mouse.init = false;
    }

    void printMat4(const glm::mat4& m4){
        std::cout<<
        m4[0][0]<<" "<<m4[1][0]<<" "<<m4[2][0]<<" "<<m4[3][0]<<"\n"<<
        m4[0][1]<<" "<<m4[1][1]<<" "<<m4[2][1]<<" "<<m4[3][1]<<"\n"<<
        m4[0][2]<<" "<<m4[1][2]<<" "<<m4[2][2]<<" "<<m4[3][2]<<"\n"<<
        m4[0][3]<<" "<<m4[1][3]<<" "<<m4[2][3]<<" "<<m4[3][3]<<"\n";
    }

    void mat2str(char* str,const glm::mat4& m4 ){
        sprintf(str, "%2.2f %2.2f %2.2f %2.2f\n%2.2f %2.2f %2.2f %2.2f\n%2.2f %2.2f %2.2f %2.2f\n%2.2f %2.2f %2.2f %2.2f\n",
        m4[0][0],m4[1][0],m4[2][0],m4[3][0],
        m4[0][1],m4[1][1],m4[2][1],m4[3][1],
        m4[0][2],m4[1][2],m4[2][2],m4[3][2],
        m4[0][3],m4[1][3],m4[2][3],m4[3][3]);
    }

    void load_render(){
        auto im = sv.scene.camera->get_image();
        int w=im.res.x;
        int h=im.res.y;
        uint8_t* tmp = new uint8_t[w*h*4];
        for(int i = 0; i<h;i++){
            for(int j =0; j<w; j++){
                Vec3 p = pixel_ops::linear2srgb((*im.get_pixel(j,i))*sceneGUI.exposure);
                pixel_ops::put2char(p,&tmp[(i*w+j)*4]);
                tmp[(i*w+j)*4+3] = 255;
            }
        }


        //todo: destroy should set resources to be deleted, but only delete them after their not in use
        vkDeviceWaitIdle(renderer_vk->getDevice());
        showRenderTex.im->destroy(renderer_vk->getIV().device);
        renderer_vk->createFullTextureImage(*showRenderTex.im, tmp,w,h);
        renderTexId = ImGui_ImplVulkan_AddTexture(
            showRenderTex.im->sampler,
            showRenderTex.im->view,
            showRenderTex.im->layout
        );





        delete[] tmp;
    }

    void updateMouse(){

        if(glfwGetInputMode(renderer_vk->getWindow(),GLFW_CURSOR) == GLFW_CURSOR_DISABLED){
            glfwGetCursorPos(renderer_vk->getWindow(), &mouse.x, &mouse.y);
            if(!mouse.init){
                mouse.init = true;
                mouse.old_x = mouse.x;
                mouse.old_y = mouse.y;
            }

            float mdx = mouse.x-mouse.old_x;
            float mdy = mouse.y-mouse.old_y;
            mouse.old_x = mouse.x;
            mouse.old_y = mouse.y;

            camera.controls.rotate_up(-mdy*camera_gui.mouseSpeed*0.01);
            camera.controls.rotate_right(mdx*camera_gui.mouseSpeed*0.01);
        }
    }

    void updateCamera(float dt){

        camera.controls.move_forward(dt*camera_gui.movementSpeed*camera_gui.move.z);
        camera.controls.move_right(dt*camera_gui.movementSpeed*camera_gui.move.x);
        camera.controls.move_up(dt*camera_gui.movementSpeed*camera_gui.move.y);
        camera.controls.tilt_right(dt*camera_gui.tiltSpeed*camera_gui.rotate.z);
        camera.fov = 2.0*atan(camera_gui.tan_fov);
    }


    //sync rules
    /*
    - Only one other thread may have access to local variables.
    - Before starting the thread, set current_task != NONE
    - In the thread, if not exiting, NEVER set current_task = NONE
    - If current_task = NONE, then it can never asynchronously be set to
        != NONE

    shared vars need to check if current_task == NONE before accessing
    NONE is the only safe task name, as all others can be changed asynchronously.
    */


    void frame(){

        static char display_str[2*MAX_CHARS] = "";
        if(!clock.isRunning()){
            clock.start();
        }

        clock.tick();
        float dt = 0.016;

        updateMouse();
        updateCamera(dt);
;
        auto m = camera.get_pose();
        mat2str(display_str,m);
        renderer_vk->getCamera() = {m, camera.fov};


        unqueueAndRunTasks();

        //if no tasks are running
        if(reloadScene && current_task==Task::Name::NONE){
            loadSceneToRenderer();
            reloadScene = false;
        }

        if(tasks.size()==0){
            if(sv.load_render){
                load_render();
            }

            sv.load_render = false;
        }

        ImGui::Begin("Render Result");

        float aspect = float(showRenderTex.im->width)/float(showRenderTex.im->height);
        float height = ImGui::GetWindowContentRegionMax().y-ImGui::GetWindowContentRegionMin().y;

        ImGui::Image((ImTextureID)renderTexId,
         {height*aspect,
         height
         });
        ImGui::End();



        ImGui::Begin("CPP-PT, a pathtracer in cpp");

        ImGui::Text(display_str);

        //if we change the input
        if(ImGui::InputText("Scene Filename",sceneGUI.scene_filename,MAX_CHARS)){
            sceneUpToDate = false;
        }

        loadSceneButton();


        ImGui::InputText("Output Name",sceneGUI.output_name,MAX_CHARS);
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
        camera_gui.imgui_panel();

        ImGui::DragInt("Res X", &sceneGUI.res.x);
        ImGui::DragInt("Res Y", &sceneGUI.res.y);
        ImGui::DragInt("SPP", &sceneGUI.spp);
        ImGui::DragFloat("Viewport Exposure", &sceneGUI.exposure,0.01,10.0);

        if(ImGui::Button("Render") && current_task == Task::Name::NONE){
            render_scene();
            std::cout<<"not blocked"<<std::endl;
        }


        showStatus();

        ImGui::End();

    }
private:

    void loadSceneButton(){

        bool local_help = sceneUpToDate;

        if(local_help){
            ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
        }
        if(ImGui::Button("Load Scene") && current_task == Task::Name::NONE && !local_help){

            sv.rs.scene_name = std::string(sceneGUI.scene_filename);
            tasks.push({

                [&](){
                    load();
                    sceneUpToDate = true;
                    reloadScene = true;
                },
                Task::Name::LOADING
            });

        }
        if(local_help){
            ImGui::PopStyleVar();
        }
    }
    void unqueueAndRunTasks(){
        if(tasks.size()>0 && current_task == Task::Name::NONE){
            if(task_thread!=nullptr){
                tasks.pop();
                std::cout<<"Joining"<<std::endl;
                task_thread->join();
                task_thread=nullptr;
            }
            if(tasks.size()>0){
                current_task = tasks.front().name;
                task_thread = std::make_unique<std::thread>([&]{
                    tasks.front().task();
                    current_task = Task::Name::NONE;
                });
                //TODO: when the app closes witha  hanging thread, what happens?
            }

        }
    }

    void loadSceneToRenderer(){
        VickyR::Scene s;

        for(auto m: sv.scene_data.meshes){
            std::vector<uint32_t> idx;
            std::vector<VickyR::Vertex> verts;
            for(auto t: m->getTriangles()){
                idx.push_back(t.at(0));
                idx.push_back(t.at(1));
                idx.push_back(t.at(2));
            }
            for(int v = 0; v<m->n_vertices(); v++){

                Vec3 p = m->get_vertex(v);
                Vec3 n= m->get_normal(v);
                Vec2 uv = m->get_uv(v);

                verts.push_back({
                   glm::vec3(p.x,p.y,p.z),
                   glm::vec3(n.x,n.y,n.z),
                   glm::vec2(uv.x,uv.y)
                });
            }
            s.addObject(std::move(verts),std::move(idx));
        }


        renderer_vk->getScene() = s;

        try{
            //TODO: waht about no scene?
            auto c = std::dynamic_pointer_cast<CameraPerspective>(sv.scene.camera);

            Vec3 r = c->getCoords().col(0);
            Vec3 b = c->getCoords().col(2);



            glm::vec3 backward = glm::vec3(
                                    b.x,
                                    b.y,
                                    b.z
                                );

            glm::vec3 right =  glm::vec3(
                                    r.x,
                                    r.y,
                                    r.z
                                );

            glm::vec3 o = glm::vec3(
                                    c->getOrigin().x,
                                    c->getOrigin().y,
                                    c->getOrigin().z
                                );



            camera = DebugTools::Camera<DebugTools::WalkControls>(
                DebugTools::WalkControls::from_vecs(
                    -backward,
                    right,
                    glm::vec3(0.0f,0.0f,1.0f),
                    o
                ),
                c->getFovy()
            );

            sceneGUI.res.x = c->get_image().res.x;
            sceneGUI.res.y = c->get_image().res.y;

            camera_gui.tan_fov = tan(0.5*c->getFovy());

        } catch (std::bad_cast e){
            std::cout<<e.what()<<std::endl;
        }


        renderer_vk->loadScene();
    }

    void showStatus(){
        if(current_task == Task::Name::LOADING){
            ImGui::Text("Loading Scene");
        }

        if(current_task == Task::Name::RENDERING && sv.renderer){
            //race condition?!!?
            ImGui::Text("Rendering Scene");


            sv.renderer->updateProgress(status);

            ImGui::ProgressBar(status.percentage);
            ImGui::SameLine(0.0f, ImGui::GetStyle().ItemInnerSpacing.x);
            char buf[32];
            sprintf(buf, "Task %02d/%02d",status.current_task, status.total_tasks);
            //std::cout<<buf<<" "<<status.percentage<<" "<<status.task_description<<std::endl;

            ImGui::Text(buf);
            ImGui::Text(status.task_description);
        }

    }

    void sync(){
        sv.rs.resX = sceneGUI.res.x;
        sv.rs.resY = sceneGUI.res.y;
        sv.rs.renderer = current_renderer;
        sv.rs.spp = sceneGUI.spp;
        sv.rs.scene_name = std::string(sceneGUI.scene_filename);
        sv.rs.output_name = sceneGUI.output_name;
        glm::vec3 o = camera.controls.origin;
        glm::vec3 f = camera.controls.forward;
        glm::vec3 u = glm::normalize(glm::cross(camera.controls.right,camera.controls.forward));
        sv.scene.camera = std::make_shared<CameraPerspective>(
            Vec2i(sv.rs.resX,sv.rs.resY),
            float(tan(0.5f*camera.fov)),
            Vec3(o.x,o.y,o.z),
            Vec3(f.x,f.y,f.z),
            Vec3(u.x,u.y,u.z)
        );
    }
    void render_scene(){
        sync();
        sv.renderer=reg->create(current_renderer,sv.rs);

        tasks.push({
            [&](){

            std::string renderer_name = reg->names[current_renderer];

            std::stringstream ss;
            ss<<sv.rs.output_name<<"x"<<sv.rs.resX<<"y"
            <<sv.rs.resY<<"spp"<<sv.rs.spp
            <<renderer_name;



            //global var access
            sv.rs.output_name = ss.str();

            //global var access
            if(!sceneUpToDate){
                load();
            }

            //global var access
            current_task = Task::Name::RENDERING;

            //global var access
            sv.renderer->setRendererMessenger(&messenger);
            sv.renderer->renderW(sv.scene, sv.rs.output_name);
            sv.load_render = true;
            },
            Task::Name::LOADING
        });
    }
    void load(){
        sync();
        sv.scene = Scene();
        sv.scene_data = SceneData();
        Loader::load_scene(&sv.scene,&sv.scene_data,&sv.rs,sv.rs.scene_name);


    }
};

void RenderMessengerGUI::eval() const  {
    rg->current_task = Task::Name::NONE;
}


}
}
#endif