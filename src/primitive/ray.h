#ifndef CPPPT_RAY
#define CPPPT_RAY


#include <math/math.h>
#include <limits>
//#define RAY_STATISTICS

namespace cpppt{
struct Ray{
    
    Vec3 o;
    Vec3 d;
    float max_t;

    #ifdef RAY_STATISTICS
    int tests = 0;
    #endif
    
    Ray(const Vec3& o=Vec3(0.0,0.0,0.0), const Vec3& d = Vec3(1.0,0.0,0.0)): o(o),d(d), max_t(std::numeric_limits<float>::infinity()) {}

};

}

#endif