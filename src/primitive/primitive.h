#ifndef CPPPT_PRIMITIVE
#define CPPPT_PRIMITIVE

#include <primitive/ray.h>
#include <shape/intersection.h>
namespace cpppt{


class Primitive {
public:
    virtual bool intersect(Ray& r, Intersection* is) const  = 0;
    virtual bool intersectAny(Ray& r) const  = 0;
};

}
#endif