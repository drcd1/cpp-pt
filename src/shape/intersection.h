#ifndef CPPPT_INTERSECTION
#define CPPPT_INTERSECTION

#include <math/math.h>
#include <bxdf/material.h>

namespace cpppt{

struct Intersection{
    Vec3 texture_coords;
    Vec3 bitangent;
    Vec3 tangent;
    Vec3 normal;
    Vec3 hitpoint;
    std::shared_ptr<Material> material;
};


}

#endif