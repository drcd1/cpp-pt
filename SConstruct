
VULKAN_DIR = r"D:\cpp_libs\VulkanSDK\1.2.170.0"
GLFW_DIR = r"D:\cpp_libs\glfw-3.3.4.bin.WIN64"
GLM_DIR = r"D:\cpp_libs\glm-0.9.9.8\glm"
STB_DIR = r"D:\cpp_libs\stb"
GLFW_LIB_DIR = "lib-vc2019"
IMGUI_DIR = r"D:\cpp_libs\imgui"
VICKYR_DIR = r"D:\VickyR"

import os

include_dir = [os.path.join(VULKAN_DIR,"include"),
               os.path.join(GLFW_DIR,"include"),
               GLM_DIR,
               STB_DIR,
               IMGUI_DIR,
               os.path.join(VICKYR_DIR,"src"),
               os.path.join(IMGUI_DIR, "backends")]

lib_dir = [os.path.join(VULKAN_DIR, "lib"), os.path.join(GLFW_DIR,GLFW_LIB_DIR)]
libs = ["vulkan-1","glfw3"]

windows_libs = [
"kernel32.lib",
"user32.lib",
"gdi32.lib",
"winspool.lib",
"comdlg32.lib",
"advapi32.lib",
"shell32.lib",
"ole32.lib",
"oleaut32.lib",
"uuid.lib",
"odbc32.lib",
"odbccp32.lib"
]




env = Environment(CPPPATH = ['src'] + include_dir)
debug = ARGUMENTS.get('debug', 0)
no_gui = ARGUMENTS.get('no_gui', 0)

if (int(debug)):
    env.Append(CCFLAGS = '/Od /openmp /std:c++17 /MD /EHsc')
    env.Append(PDB = 'build/main.pdb')
    env.Append(CCPDBFLAGS = '/Zi')
    env.Append(LINKFLAGS = '/DEBUG')
else:
    env.Append(CCFLAGS = "/O2 /openmp /std:c++17 /MD /EHsc")



source_files = ["src/main.cpp","src/debug_tools/app.cpp"]


env.Append(LIBPATH=lib_dir)
env.Append(LIBS = libs + windows_libs)

lib_source_files = [os.path.join(IMGUI_DIR,"backends/imgui_impl_glfw.cpp"),
                    os.path.join(IMGUI_DIR,"backends/imgui_impl_vulkan.cpp"),
                    os.path.join(IMGUI_DIR,"imgui.cpp"),
                    os.path.join(IMGUI_DIR,"imgui_draw.cpp"),
                    os.path.join(IMGUI_DIR,"imgui_demo.cpp"),
                    os.path.join(IMGUI_DIR,"imgui_tables.cpp"),
                    os.path.join(IMGUI_DIR,"imgui_widgets.cpp")]
if(not int(no_gui)):
    env.Append(CPPDEFINES=['USE_GUI'])

env.Program("build/renderer",source_files + lib_source_files)





#### Build shaders
#shader_builder = Builder(action = 'D:/cpp_libs/VulkanSDK/1.2.170.0/Bin32/glslc.exe $SOURCE -o $TARGET')
#env = Environment()
#env.Append(BUILDERS={'Shader' : bld})
#for f in os.listdir("shaders"):
#    fn = os.path.join("shaders",f)
#    if os.path.isfile(fn):
#        env.Shader( os.path.join("build/shaders",f+".spv"), fn)