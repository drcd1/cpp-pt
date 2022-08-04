#ifndef CPPPT_PRIMITIVE_LEAF
#define CPPPT_PRIMITIVE_LEAF

#include <primitive/ray.h>
#include <primitive/primitive_leaf.h>
#include <primitive/aabb.h>
#include <shape/shape.h>
#include <memory.h>
namespace cpppt{


class PrimitiveLeaf: public Primitive {
private:
    std::shared_ptr<Shape> shape;
    std::shared_ptr<BxDF> bxdf;
    AABB aabb;
public:
    PrimitiveLeaf(std::shared_ptr<Shape> shape, std::shared_ptr<BxDF> bxdf): 
    shape(shape), bxdf(bxdf),
    aabb(shape->get_bounds()){}
    bool intersect(Ray& r, Intersection* is) const {
        bool intersected = shape->intersect(r,is);
        if(intersected){
            is->material = bxdf;
        }
        return intersected;
    }
    bool intersectAny(Ray& r) const {
        return shape->intersectAny(r);
    }
    AABB get_bounds() const {
        return aabb;
    }
};


}

#endif