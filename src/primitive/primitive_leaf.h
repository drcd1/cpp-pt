#ifndef CPPPT_PRIMITIVE_LEAF
#define CPPPT_PRIMITIVE_LEAF

#include <primitive/ray.h>
#include <primitive/primitive_leaf.h>
#include <primitive/aabb.h>
#include <shape/shape.h>
#include <memory.h>
#include <bxdf/standard_material.h>


namespace cpppt{

class ShapeLight;

//TODO. make sure primitive lifetime is the same as lights
class PrimitiveLeaf: public Primitive {
private:
    std::shared_ptr<Shape> shape;
    std::shared_ptr<Material> material;
    ShapeLight* light;
    AABB aabb;
public:
    PrimitiveLeaf(std::shared_ptr<Shape> shape, std::shared_ptr<Material> material):
    shape(shape), material(material),light(nullptr),
    aabb(shape->get_bounds()){}
    bool intersect(Ray& r, Intersection* is) const {
        bool intersected = shape->intersect(r,is);
        if(intersected){
            is->material = material;
            is->primitive = this;
        }
        return intersected;
    }
    bool intersect_any(Ray& r) const {
        return shape->intersect_any(r);
    }
    AABB get_bounds() const {
        return aabb;
    }

    void set_light(ShapeLight* sl) {
        light = sl;
    }

    const ShapeLight* get_light() const {
        return light;
    }

    const std::shared_ptr<Shape> get_shape() const {
        return shape;
    }

    const std::shared_ptr<Material> get_material() const {
        return material;
    }
};


}

#endif