#ifndef CPPPT_SIMPLEGROUP
#define CPPPT_SIMPLEGROUP

#include <primitive/ray.h>
#include <memory>


namespace cpppt{


class SimpleGroup: public Primitive {
private:
    std::vector<std::shared_ptr<Primitive>> primitives;
    AABB aabb;

public:    

    void add(std::shared_ptr<Primitive> primitive){
        primitives.push_back(primitive);
    }

    bool intersect(Ray& r, Intersection* is) const {
        bool intersected = false;
        for(int i = 0; i<primitives.size(); i++){
            bool tmp = primitives.at(i)->intersect(r,is);
            intersected =  tmp || intersected;
        }
        return intersected;
    }
    bool intersectAny(Ray& r) const {
        bool intersected = false;
        for(int i = 0; i<primitives.size(); i++){
            intersected = primitives.at(i)->intersectAny(r);
            if(intersected){
                return intersected;
            }
        }
        return false;
    }

    AABB get_bounds() const {
        return aabb;
    }
};


}

#endif