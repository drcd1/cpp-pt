        
#ifndef CPPPT_TRIANGLE
#define CPPPT_TRIANGLE


namespace cpppt {   
    class Triangle: public Shape {
        public:
        virtual bool intersect(Ray& r, Intersection* it) const = 0;
        virtual bool intersectAny(Ray& r) const = 0;

    };
}

#endif