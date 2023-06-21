#ifndef CPPPT_DEBUGTOOLS_TYPES
#define CPPPT_DEBUGTOOLS_TYPES

#include <math/math.h>
#include <vector>
namespace cpppt {

    template <typename T>
    using GUIAnimation = std::vector<T>;  
    
    typedef std::vector<Vec3> GUILightpath;
}

#endif