#ifndef CPPPT_DEBUG_TOOLS_IMAGE
#define CPPPT_DEBUG_TOOLS_IMAGE

#include <image/rgb_image.h>
#include <vickyr.h>
#include <imgui.h>
#include <imgui_impl_vulkan.h>
namespace cpppt{
namespace DebugTools{

struct DisplayTexture{
    RgbImage image;
    VickyR::VickyR* renderer;
    VickyR::AllocatedTextureImage tex;
    VkDescriptorSet renderTexId;

    DisplayTexture(RgbImage& a_image, VickyR::VickyR* a_renderer):image(a_image),renderer(a_renderer){
        int w=image.res.x;
        int h=image.res.y;
        uint8_t* tmp = new uint8_t[w*h*4];
        for(int i = 0; i<h;i++){
            for(int j =0; j<w; j++){
                Vec3 p = pixel_ops::linear2srgb((*image.get_pixel(j,i)));
                pixel_ops::put2char(p,&tmp[(i*w+j)*4]);
                tmp[(i*w+j)*4+3] = 255;
            }
        }
        tex =  renderer->getIV().imagePool.allocate();



        //showRenderTex.im->destroy(renderer->getIV().device);
        renderer->createFullTextureImage(*tex.im, tmp,w,h);
        renderTexId = ImGui_ImplVulkan_AddTexture(
            tex.im->sampler,
            tex.im->view,
            tex.im->layout
        );
    }

    ~DisplayTexture(){
         //todo: destroy should set resources to be deleted, but only delete them after their not in use
        vkDeviceWaitIdle(renderer->getDevice());
        tex.im->destroy(renderer->getDevice());
        renderer->getIV().imagePool.deallocate(tex.id);
    }

}

class ImageGUI{
    std::vector
}


}
}


#endif