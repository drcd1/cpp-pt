#ifndef CPPPT_PRIMITIVE
#define CPPPT_PRIMITIVE

#include <primitive/ray.h>
#include <shape/intersection.h>
#include <primitive/aabb.h>
namespace cpppt{


class Primitive {
public:
    virtual bool intersect(Ray& r, Intersection* is) const  = 0;
    virtual bool intersect_any(Ray& r) const  = 0;
    virtual AABB get_bounds() const = 0;
};



}
#endif