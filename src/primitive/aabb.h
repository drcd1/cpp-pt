#ifndef CPPPT_AABB
#define CPPPT_AABB

#include <shape/intersection.h>
#include <primitive/ray.h>


namespace cpppt{

    class AABB{
        public:
        Vec3 min;
        Vec3 max;
        AABB(): 
            min(std::numeric_limits<float>::infinity()),
            max(-std::numeric_limits<float>::infinity())
        {}
        AABB(const Vec3& min,const Vec3& max):min(min),max(max){}
        bool intersectAny(const Ray& r) const{
            Vec3 nMin = min;
            Vec3 nMax = max;
            

            Vec3 dmin = nMin - r.o;
            Vec3 dmax = nMax - r.o;

            nMin.x = dmin.x/r.d.x;
            nMin.y = dmin.y/r.d.y;
            nMin.z = dmin.z/r.d.z;

            nMax.x = dmax.x/r.d.x;
            nMax.y = dmax.y/r.d.y;
            nMax.z = dmax.z/r.d.z;

            order2(nMin.x,nMax.x);
            order2(nMin.y,nMax.y);
            order2(nMin.z,nMax.z);

            float tmin = nMin.x;
            if(tmin<nMin.y)
                tmin = nMin.y;
            if(tmin<nMin.z && nMin.z>nMin.y)
                tmin = nMin.z;
            
            float tmax = nMax.x;
            if(tmax>nMax.y)
                tmax = nMax.y;
            if(tmax>nMax.z && nMax.z<nMax.y)
                tmax = nMax.z;

            if(tmax<tmin)
                return false;

            if(tmin>0 && tmin<r.max_t+EPS){
                return true;
            }
            if(tmax>0 /*&& tmax<r.max_t+EPS*/){
                return true;
            }
            return false;
        }
        Vec3 center() const {
            return (min+max)*0.5;
        }

        void unite(const AABB& other) {
            min = minF(min,other.min);
            max = maxF(max,other.max);
        }

        
        void unite(const Vec3& other) {
            min = minF(min,other);
            max = maxF(max,other);
        }

        float area() const {
            return 2.0*( fabs(max.x-min.x)*(max.y-min.y) +  
            fabs(max.z-min.z)*(max.y-min.y) 
            + fabs(max.x-min.x)*(max.z-min.z)
            );
        }
    };
}

#endif