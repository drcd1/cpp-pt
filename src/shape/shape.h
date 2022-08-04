#ifndef CPPPT_SHAPE
#define CPPPT_SHAPE

#include <shape/intersection.h>
#include <primitive/ray.h>
#include <primitive/aabb.h>
#include <iostream>
namespace cpppt{

    class Shape {
        public:
        virtual bool intersect(Ray& r, Intersection* it) const = 0;
        virtual bool intersectAny(Ray& r) const = 0;
        virtual AABB get_bounds() const = 0;

    };



    class Sphere: public Shape {
    private:
        Vec3 origin;
        float radius;
        public:
        Sphere(Vec3 o, float r): origin(o),radius(r){}

        AABB get_bounds() const {
            return AABB(origin-Vec3(radius),origin+Vec3(radius));
        }

        bool intersect(Ray& r, Intersection* it) const {

           Vec3 ro = r.o - origin;

            //  dot(x,x) = radius^2
            //      x = ro + r.d*t
            //  dot(x,x) = ro^2 + 2*r.d*t*ro + r.d^2*t^2
            // ro^2 - radius^2 + 2*r.d*t*ro + r.d^2*t^2 = 0
            //  a = r.d^2,   b = 2*r.d*t*ro,  c = ro^2 - radius^2

            float t0,t1;
            float a = lensqr(r.d);
            float b = 2.0*dot(r.d, ro);
            float c = lensqr(ro) - radius*radius;
            bool solved = solve_quadratic(a,b,c,&t0,&t1);
            if(!solved){
                return false;
            }
            float t;
            if(t0>EPS){
                t = t0;            
            } else if(t1>EPS){
                t = t1;
            } else {
                return false;
            }

            if(t>=r.max_t){
                return false;
            }

            r.max_t = t;

            it->hitpoint = r.o + r.d*t;
            it->normal = normalized(it->hitpoint - origin);
            it->hitpoint = it->hitpoint+ it->normal*EPS;

            //Vec3 help = it->hitpoint - origin;
            //float theta = acos(help.z/len(Vec3(help.x,help.y,0.0)));
            //float phi = atan2(help.y,help.x);

            it->texture_coords = Vec3(1.0,0.0,0.0);
            return true;
            


        }
        bool intersectAny(Ray& r) const {
            return false;
        }

    };
}
#endif