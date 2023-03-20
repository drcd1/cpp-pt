#ifndef CPPPT_DEBUG_TOOLS_CAMERA
#define CPPPT_DEBUG_TOOLS_CAMERA

#define GLM_FORCE_RADIANS
#define GLM_FORCE_DEPTH_ZERO_TO_ONE
#include <glm/vec4.hpp>
#include <glm/mat4x4.hpp>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

namespace cpppt{
namespace DebugTools{

struct CameraGUI{
    Vec3 move = Vec3(0.0,0.0,0.0);
    Vec3 rotate=Vec3(0.0,0.0,0.0);
    float tan_fov = 1.0;
    float movementSpeed=0.5;
    float mouseSpeed=0.2;
    float tiltSpeed=0.5;

    void imgui_panel(){

        ImGui::DragFloat("Movement speed", &movementSpeed,0.2, 0.0, 5.0);
        ImGui::DragFloat("Tilt speed", &tiltSpeed,0.2, 0.0, 5.0);
        ImGui::DragFloat("Mouse speed", &mouseSpeed,0.2, 0.0, 5.0);

        ImGui::DragFloat("Camera FOV", &tan_fov,0.05,0.01,10.0);
    }
};


struct WalkControls{
private:
    WalkControls(const glm::vec3 f,
    const glm::vec3 r,
    const glm::vec3 u,
    const glm::vec3 o): forward(f),right(r),up(u),origin(o){}
public:
    glm::vec3 forward;
    glm::vec3 right;
    glm::vec3 up;
    glm::vec3 origin;


    static inline WalkControls from_vecs(
        glm::vec3 f,
        glm::vec3 r,
        glm::vec3 u,
        glm::vec3 o
    ){
        return WalkControls(f,r,u,o);
    }

    static inline WalkControls from_transform(
        const glm::mat4& m,
        const glm::vec3& up
    ){
        //transform makes space into camera
        glm::vec3 forward = -glm::vec3(m[0][2],m[1][2],m[2][2]);
        glm::vec3 right =    glm::vec3(m[0][0],m[1][0],m[2][0]);
        glm::vec3 help_up = cross(right,forward);
        glm::vec3 origin = glm::vec3(m[3][0],m[3][1],m[3][2]);
        origin = origin.x*right + origin.y*help_up-origin.z*forward;

        return WalkControls(forward,right,up,origin);
    }

    void move_forward(float t){
        glm::vec3 dir = normalize(forward-dot(forward,up)*up);
        origin+=dir*t;
        //origin.z+=t;
    }

    void move_right(float t){

        glm::vec3 f_d = normalize(forward-dot(forward,up)*up);
        glm::vec3 r_d = cross(f_d,up);
        origin+=r_d*t;
        //origin.x+=t;
    }

    void move_up(float t){
        origin+=up*t;
        //origin.y+=t;
    }


    void rotate_up(float t){
        glm::vec3 f_d = normalize(forward-dot(forward,up)*up);

        glm::vec3 r_d = cross(f_d,up);
        glm::mat4 rotate = glm::rotate(glm::mat4(1.0), t, r_d);
        glm::vec3 help = (rotate*glm::vec4(forward,0.0f));
        if(dot(help, f_d)<0.01){
            return;
        }
        forward = help;
        right = (rotate*glm::vec4(right,0.0f));
    }


    void rotate_right(float t){
        glm::vec3 f_d = normalize(forward-dot(forward,up)*up);
        glm::mat4 rotate = glm::rotate(glm::mat4(1.0), -t, up);
        forward = (rotate*glm::vec4(forward,0.0f));;
        right = (rotate*glm::vec4(right,0.0f));
    }

    void tilt_right(float t){

        glm::mat4 rotate = glm::rotate(glm::mat4(1.0), t, forward);
        right = rotate*glm::vec4(right,0.0f);
    }

    glm::mat4 get_pose(){
        glm::vec3 n_up = glm::cross(right, forward);
        glm::mat4 pose = glm::mat4(
             glm::vec4(right,0.0),
             glm::vec4(n_up,0.0),
             glm::vec4(-forward,0.0),
             glm::vec4(0.0,0.0,0.0,1.0)
        );

        return glm::transpose(pose)*glm::translate(glm::mat4(1.0), -glm::vec3(origin));

       //return glm::translate(glm::mat4(1.0), -glm::vec3(origin));
    }
};

template
<typename Controls>
struct Camera{
    Controls controls;
    Camera(Controls&& controls, float fov):controls(controls),fov(fov){};

    float fov;


    glm::mat4 get_pose(){
        return controls.get_pose();
    }


};
}
}
#endif