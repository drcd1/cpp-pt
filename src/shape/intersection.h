#ifndef CPPPT_INTERSECTION
#define CPPPT_INTERSECTION

#include <math/math.h>
#include <bxdf/material.h>
#include <memory>


namespace cpppt{

class Shape;
class PrimitiveLeaf;
class BxDF;

struct Intersection{
    Vec3 texture_coords;
    Vec3 bitangent;
    Vec3 tangent;
    Vec3 normal;
    Vec3 g_normal;
    Vec3 hitpoint;
    bool computed_normal = false;
    const PrimitiveLeaf* primitive;

    std::shared_ptr<Material> material;
    std::shared_ptr<BxDF> get_bxdf(){
        if(!computed_normal){
            normal = material->get_normal_mapped(texture_coords,tangent,bitangent,normal);
            bitangent = normalized(cross(normal,tangent));
            tangent = normalized(cross(bitangent,normal));
            computed_normal = true;
        }
        return material->get_bxdf(texture_coords);
    }

};


}

#endif