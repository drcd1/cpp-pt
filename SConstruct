from config import *
import os
import re


DEBUG_TARGET = "debug"
RELEASE_TARGET = "release"
TOOLS_TARGET = "tools"
MIN_TARGET = "min"
ALL_TARGET = "all"


target = ARGUMENTS.get('target', ALL_TARGET)
debug = ARGUMENTS.get('debug', ALL_TARGET)
platform = ARGUMENTS.get('OS', Platform())

if(target == ALL_TARGET):
    targets = [MIN_TARGET, TOOLS_TARGET]
else:
    targets = [target]

if(debug==ALL_TARGET):
    debugs = [DEBUG_TARGET, RELEASE_TARGET]
else:
    debugs = [debug]


include_dir = ["src"]
lib_dir = []

include_dir_tools = []
include_dir_min = []


lib_dir_tools = []
lib_dir_min = []


if(VULKAN_DIR!=""):
    include_dir_tools += [os.path.join(VULKAN_DIR,"include")]
    lib_dir_tools +=[os.path.join(VULKAN_DIR,"lib")]

if(EMBREE_DIR!=""):
    include_dir += [os.path.join(EMBREE_DIR,"include")]
    lib_dir+=[os.path.join(EMBREE_DIR,"lib")]


if(GLFW_DIR!=""):
    include_dir_tools += [os.path.join(GLFW_DIR,"include")]
    lib_dir_tools += [os.path.join(GLFW_DIR,GLFW_LIB_DIR)]

if(GLM_DIR!=""):
    include_dir+=[GLM_DIR]
if(STB_DIR!=""):
    include_dir+=[STB_DIR]

if(VICKYR_DIR!=""):
    include_dir_tools += [os.path.join(VICKYR_DIR,"src")]
if(IMGUI_DIR!=""):
    include_dir_tools+=[IMGUI_DIR, os.path.join(IMGUI_DIR, "backends")]





libs = []
libs_tools = []
libs_min = []


def createProgram(t,d):

    var_src = "build/"+t+"_"+d+"/src/"
    var_imgui = "build/"+t+"_"+d+"/imgui/"

    VariantDir(var_src,"src",duplicate=0)
    VariantDir(var_imgui,IMGUI_DIR,duplicate=0)
    #VariantDir("build_imgui"+t+"_"+d,IMGUI_DIR)
    env = Environment(Platform=platform)


    if(t==TOOLS_TARGET):
        inc_dir = include_dir+include_dir_tools
        l_dir = lib_dir+lib_dir_tools
        ls = libs + libs_tools
    elif(t==MIN_TARGET):
        inc_dir = include_dir+include_dir_min
        l_dir = lib_dir+lib_dir_min
        ls = libs + libs_min

    env.Append(CPPPATH=inc_dir)
    env.Append(LIBPATH=l_dir)
    env.Append(LIBS = ls)


    if(platform.name=="win32"):
        env.Append(CCFLAGS = ['/EHsc','/openmp', '/MD' ,'/std:c++17'])
    else : #linux :
        env.Append(LINKFLAGS=["-fopenmp"])
        env.Append(CCFLAGS= ['-std=c++17','-fopenmp' ])

    if(platform.name=="win32"):
        if d==DEBUG_TARGET:
            print("debug build!");
            env.Append(PDB = 'debug/renderer_'+t+"_"+d+'.pdb')
            env.Append(CCPDBFLAGS = '/Zi')
        else:
            env.Append(CCFLAGS = '/O2')
    else :#linux:
        if d==DEBUG_TARGET:
            print("debug build!");
            env.Append(CCFLAGS='-g')
        else:
            env.Append(CCFLAGS = '-O2')

    env.Append(LIBPATH=lib_dir)
    env.Append(LIBS = libs)

    lib_source_files = []
    main_file = ""
    source_files = []
    if(t==TOOLS_TARGET):
        main_file=var_src+'debug_tools/app.cpp'

    elif(t== MIN_TARGET):
        main_file= var_src+'main.cpp'

    if(t == TOOLS_TARGET):
        lib_source_files += [os.path.join(var_imgui,"backends/imgui_impl_glfw.cpp"),
                        os.path.join(var_imgui,"backends/imgui_impl_vulkan.cpp"),
                        os.path.join(var_imgui,"imgui.cpp"),
                        os.path.join(var_imgui,"imgui_draw.cpp"),
                        os.path.join(var_imgui,"imgui_demo.cpp"),
                        os.path.join(var_imgui,"imgui_tables.cpp"),
                        os.path.join(var_imgui,"imgui_widgets.cpp")]

        vd_def = ("VICKYR_DIR=\\\""+re.escape(re.escape(VICKYR_DIR+"build/"))+"\\\"")
        env.Append(CPPDEFINES=['USE_GUI',vd_def])
        #source_files += ["src/debug_tools/app.cpp"]

    files = [main_file]+source_files+lib_source_files
    print("asda"+t+d)
    print(files)
    env.Program("build/renderer_"+t+"_"+d,files)

print("PLatform is ",platform)

if(platform.name=="win32"):
    libs_tools += ["vulkan-1","glfw3","gdi32"]

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
    "odbccp32.lib",
    "embree4.lib",
    "tbb.lib"
    ]
    libs_tools += windows_libs
    libs_min += ["embree4.lib","tbb.lib"]
else: #linux
    libs_tools += ["glfw","vulkan","dl","X11","Xxf86vm","Xrandr","Xi"]
    libs += ["pthread","embree4","tbb"]



for t in targets:
    for d in debugs:
        print(t,d)
        createProgram(t,d)
