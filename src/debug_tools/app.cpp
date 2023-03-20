

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_vulkan.h"
#include <debug_tools/gui/renderer.h>
#include <vickyr.h>

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image.h>
#include <stb_image_write.h>

cpppt::DebugTools::RendererGUI* g_r_ptr =nullptr;

void mouseButtonCallback(GLFWwindow* window,int button, int action, int mods){
	if(ImGui::GetIO().WantCaptureMouse || g_r_ptr == nullptr){
		return;
	}
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS){
		glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
		g_r_ptr->resetMouse();
		if (glfwRawMouseMotionSupported())
    		glfwSetInputMode(window, GLFW_RAW_MOUSE_MOTION, GLFW_TRUE);
	}
}

void keyCallback(GLFWwindow* window,int button, int scancode, int action, int mods){
	if(ImGui::GetIO().WantCaptureMouse){
		return;
	}
	if(button==GLFW_KEY_ESCAPE && action == GLFW_PRESS){
		if(glfwGetInputMode(window,GLFW_CURSOR) == GLFW_CURSOR_DISABLED){
			glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
			if (glfwRawMouseMotionSupported())
    			glfwSetInputMode(window, GLFW_RAW_MOUSE_MOTION, GLFW_FALSE);
		}
	}

	if(g_r_ptr){
	if(action==GLFW_PRESS){
		if(button==GLFW_KEY_W){
			g_r_ptr->getCameraGUI().move.z = 1.0;
		}
		if(button==GLFW_KEY_S){
			g_r_ptr->getCameraGUI().move.z = -1.0;
		}

		if(button==GLFW_KEY_A){
			g_r_ptr->getCameraGUI().move.x = -1.0;
		}

		if(button==GLFW_KEY_D){
			g_r_ptr->getCameraGUI().move.x = 1.0;
		}

		if(button==GLFW_KEY_Z){
			g_r_ptr->getCameraGUI().move.y = -1.0;
		}
		if(button==GLFW_KEY_X){
			g_r_ptr->getCameraGUI().move.y = 1.0;
		}

		if(button==GLFW_KEY_Q){
			g_r_ptr->getCameraGUI().rotate.z = -1.0;
		}
		if(button==GLFW_KEY_E){
			g_r_ptr->getCameraGUI().rotate.z = 1.0;
		}
	} else if(action==GLFW_RELEASE){
		if(button==GLFW_KEY_W){
			g_r_ptr->getCameraGUI().move.z = 0.0;
		}
		if(button==GLFW_KEY_S){
			g_r_ptr->getCameraGUI().move.z = 0.0;
		}

		if(button==GLFW_KEY_A){
			g_r_ptr->getCameraGUI().move.x = 0.0;
		}

		if(button==GLFW_KEY_D){
			g_r_ptr->getCameraGUI().move.x = 0.0;
		}

		if(button==GLFW_KEY_Z){
			g_r_ptr->getCameraGUI().move.y = 0.0;
		}
		if(button==GLFW_KEY_X){
			g_r_ptr->getCameraGUI().move.y = 0.0;
		}

		if(button==GLFW_KEY_Q){
			g_r_ptr->getCameraGUI().rotate.z = 0.0;
		}
		if(button==GLFW_KEY_E){
			g_r_ptr->getCameraGUI().rotate.z = 0.0;
		}
	}
	}
}


void init_callbacks(GLFWwindow* window){
	glfwSetMouseButtonCallback(window, mouseButtonCallback);
	glfwSetKeyCallback(window, keyCallback);
	glfwSetMouseButtonCallback(window, mouseButtonCallback);
	glfwSetMouseButtonCallback(window, mouseButtonCallback);
}

void init_imgui(VickyR::VickyR& r){

	init_callbacks(r.getWindow());
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

	//todo: Hook up renderer callbacks
	//before this
	//this initializes imgui for Glfw
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
    VickyR::VickyR renderer;
    renderer.setInitHook([&](){
        init_imgui(renderer);
    });


    renderer.init();
    // Setup Dear ImGui context
    //IMGUI_CHECKVERSION();
    //ImGui::CreateContext();
    //ImGuiIO& io = ImGui::GetIO(); (void)io;

    //ImGui::StyleColorsDark();


    cpppt::DebugTools::RendererGUI rgui(renderer);

    renderer.setRenderPassHook([&](VkCommandBuffer& cmd){
        ImGui_ImplVulkan_RenderDrawData(ImGui::GetDrawData(), cmd);
    });

	g_r_ptr = &rgui;




    //run frame
    renderer.mainLoop([&](){
        ImGui_ImplVulkan_NewFrame();
		ImGui_ImplGlfw_NewFrame();

		ImGui::NewFrame();
        //imgui commands

		rgui.frame();
        ImGui::ShowDemoWindow();
        ImGui::Render();
    });
    renderer.cleanup();
}