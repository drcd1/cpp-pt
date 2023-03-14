

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_vulkan.h"
#include <debug_tools/gui/renderer.h>
#include <vickyr.h>

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image.h>
#include <stb_image_write.h>



void init_imgui(VickyR& r){
	VkDescriptorPoolSize pool_sizes[] =
	{
		{ VK_DESCRIPTOR_TYPE_SAMPLER, 1000 },
		{ VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1000 },
		{ VK_DESCRIPTOR_TYPE_SAMPLED_IMAGE, 1000 },
		{ VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, 1000 },
		{ VK_DESCRIPTOR_TYPE_UNIFORM_TEXEL_BUFFER, 1000 },
		{ VK_DESCRIPTOR_TYPE_STORAGE_TEXEL_BUFFER, 1000 },
		{ VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1000 },
		{ VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1000 },
		{ VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER_DYNAMIC, 1000 },
		{ VK_DESCRIPTOR_TYPE_STORAGE_BUFFER_DYNAMIC, 1000 },
		{ VK_DESCRIPTOR_TYPE_INPUT_ATTACHMENT, 1000 }
	};

	VkDescriptorPoolCreateInfo pool_info = {};
	pool_info.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
	pool_info.flags = VK_DESCRIPTOR_POOL_CREATE_FREE_DESCRIPTOR_SET_BIT;
	pool_info.maxSets = 1000;
	pool_info.poolSizeCount = std::size(pool_sizes);
	pool_info.pPoolSizes = pool_sizes;

	VkDescriptorPool imguiPool;
	if(vkCreateDescriptorPool(r.getDevice(), &pool_info, nullptr, &imguiPool)){
        throw std::runtime_error("Error creating descriptor pool");
    };


	// 2: initialize imgui library

	//this initializes the core structures of imgui
	ImGui::CreateContext();

	//this initializes imgui for SDL
	ImGui_ImplGlfw_InitForVulkan(r.getWindow(),true);

	//this initializes imgui for Vulkan
	ImGui_ImplVulkan_InitInfo init_info = {};
	init_info.Instance = r.getInstance();
	init_info.PhysicalDevice = r.getPhysicalDevice();
	init_info.Device = r.getDevice();
	init_info.Queue = r.getGraphicsQueue();
	init_info.DescriptorPool = imguiPool;
	init_info.MinImageCount = r.ic;
	init_info.ImageCount = r.ic;
	init_info.MSAASamples = VK_SAMPLE_COUNT_1_BIT;

	ImGui_ImplVulkan_Init(&init_info, r.getRenderPass());

	//execute a gpu command to upload imgui font textures
	VkCommandBuffer commandBuffer = r.getIV().beginSingleTimeCommands();
	ImGui_ImplVulkan_CreateFontsTexture(commandBuffer);
	r.getIV().endSingleTimeCommands(commandBuffer);

	//clear font textures from cpu data
	ImGui_ImplVulkan_DestroyFontUploadObjects();

    VkDevice d = r.getDevice();

	//add the destroy the imgui created structures
	r.addCleanup([=]() {
		vkDestroyDescriptorPool(d, imguiPool, nullptr);
		ImGui_ImplVulkan_Shutdown();
	});
}

int main(){
    VickyR renderer;
    renderer.setInitHook([&](){
        init_imgui(renderer);
    });

    renderer.init();
    // Setup Dear ImGui context
    //IMGUI_CHECKVERSION();
    //ImGui::CreateContext();
    //ImGuiIO& io = ImGui::GetIO(); (void)io;

    //ImGui::StyleColorsDark();

    renderer.setRenderPassHook([&](VkCommandBuffer& cmd){
        ImGui_ImplVulkan_RenderDrawData(ImGui::GetDrawData(), cmd);
    });



    //run frame
    renderer.mainLoop([&](){
        ImGui_ImplVulkan_NewFrame();
		ImGui_ImplGlfw_NewFrame();

		ImGui::NewFrame();
        //imgui commands
        ImGui::ShowDemoWindow();
        ImGui::Render();
    });
    renderer.cleanup();
}