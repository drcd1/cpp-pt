#ifndef CPPPT_PRIMITIVE_LEAF
#define CPPPT_PRIMITIVE_LEAF

#include <primitive/ray.h>
#include <primitive/primitive_leaf.h>
#include <shape/shape.h>
#include <memory.h>
namespace cpppt{


class PrimitiveLeaf: public Primitive {
private:
    std::shared_ptr<Shape> shape;
    std::shared_ptr<BxDF> bxdf;

public:
    PrimitiveLeaf(std::shared_ptr<Shape> shape, std::shared_ptr<BxDF> bxdf): shape(shape), bxdf(bxdf){}
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
};


}

#endif