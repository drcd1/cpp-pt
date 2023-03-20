
#ifndef CPPPT_DEBUG_TOOLS_SCENE
#define CPPPT_DEBUG_TOOLS_SCENE

#include <math/math.h>
#include <debug_tools/gui/defs.h>

namespace cpppt{
    namespace DebugTools{


struct SceneGUI{
    char scene_filename[MAX_CHARS] = "../scenes/scene_hdr_real/scene.ptscene";
    char output_name[MAX_CHARS] = "scene_hdr";
    Vec2i res={256,256};
    int spp=1;
    float exposure = 1.0;
};
}
}
#endif