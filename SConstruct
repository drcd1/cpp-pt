
#VULKAN_DIR = r"D:\cpp_libs\VulkanSDK\1.2.170.0"
#GLFW_DIR = r"D:\cpp_libs\glfw-3.3.4.bin.WIN64"
#GLM_DIR = r"D:\cpp_libs\glm-0.9.9.8\glm"
#STB_DIR = r"D:\cpp_libs\stb"
#GLFW_LIB_DIR = "lib-vc2019"
#IMGUI_DIR = r"D:\cpp_libs\imgui"
#VICKYR_DIR = r"D:\VickyR"

VULKAN_DIR = r""
GLFW_DIR = r""
GLM_DIR = r""
STB_DIR = r"/home/duarte/stb/"
IMGUI_DIR = r"/home/duarte/cpplibs/imgui"
GLFW_LIB_DIR = ""
VICKYR_DIR = r"/home/duarte/VickyR"

import os
include_dir = ['src']

include_dir = ["src"]
lib_dir = []


if(VULKAN_DIR!=""):
    include_dir += [os.path.join(VULKAN_DIR,"include")]
    lib_dir+=[os.path.join(VULKAN_DIR,"lib")]

if(GLFW_DIR!=""):
    include_dir += os.path.join(GLFW_DIR,"include")
    lib_dir += [os.path.join(GLFW_DIR,GLFW_LIB_DIR)]

if(GLM_DIR!=""):
    include_dir+=[GLM_DIR]
if(STB_DIR!=""):
    include_dir+=[STB_DIR]

if(VICKYR_DIR!=""):
    include_dir += [os.path.join(VICKYR_DIR,"src")]
if(IMGUI_DIR!=""):
    include_dir+=[IMGUI_DIR, os.path.join(IMGUI_DIR, "backends")]


debug = ARGUMENTS.get('debug', 0)
no_gui = ARGUMENTS.get('no_gui', 0)
platform = ARGUMENTS.get('OS', Platform())


if(platform=="win32"):
    libs = ["vulkan-1","glfw3","gdi32"]

#if using windows
if(platform == "win32"):
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
    libs = libs+windows_libs
else: #linux 
    libs = ["glfw","vulkan","dl","pthread","X11","Xxf86vm","Xrandr","Xi"]

env = Environment(Platfor=platform)

env.Append(CPPPATH=include_dir)
env.Append(LIBPATH=lib_dir)
env.Append(LIBS = libs)

if(platform=="win32"):
    env.Append(CCFLAGS = ['/EHsc','/openmp', '/MD' ,'/std:c++17'])
else : #linux :
    env.Append(LINKFLAGS=["-fopenmp"])
    env.Append(CCFLAGS= ['-std=c++17','-fopenmp' ])

if(platform=="win32"):
    if int(debug):
        print("debug build!");
        env.Append(PDB = 'debug/main.pdb')
        env.Append(CCPDBFLAGS = '/Zi')
    else:
        env.Append(CCFLAGS = '/O2')
else :#linux:
    if int(debug):
        print("debug build!");
        env.Append(CCFLAGS='-g')
    else:
        env.Append(CCFLAGS = '-O2')
    







source_files = ["src/main.cpp","src/debug_tools/app.cpp"]


env.Append(LIBPATH=lib_dir)
env.Append(LIBS = libs)

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