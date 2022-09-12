
#ifndef CPPPT_TRIANGLE
#define CPPPT_TRIANGLE

#include <math/math.h>
#include <shape/intersection.h>
#include <shape/shape.h>
#include <shape/mesh.h>

namespace cpppt {
    class Triangle: public Shape {
        private:
            const Mesh* mesh;
            int id;
        public:
        Triangle(const Mesh* mesh, int id): mesh(mesh),id(id){
        }
        bool intersect(Ray& r, Intersection* it) const {
            Vec3 a,b,c;
            std::array<int,3> tri =  mesh->get_triangle(id);
            a = mesh->get_vertex(tri[0]);
            b = mesh->get_vertex(tri[1]);
            c = mesh->get_vertex(tri[2]);

            //p = a + e1*x+e2*y;
            //p = o +d*t

            //o + d*t = a + e1*x+e2*y
            //o - a = -d*t + e1*x + e2*y

            /*
                       [x]
            [e1 e2 -d] [y] = [o - a]
                       [t]
            */
            Vec3 v1 = b-a;
            Vec3 v2 = c-a;
            Vec3 v3 = -r.d;
            Vec3 vs = (Vec3) (r.o - a);

            float invDenominator = 1.0/(v1.x * v2.y * v3.z +
                                    v2.x * v3.y * v1.z +
                                    v3.x * v1.y * v2.z -
                                    v3.x * v2.y * v1.z -
                                    v3.y * v2.z * v1.x -
                                    v3.z * v2.x * v1.y);

            float t = (v1.x * v2.y * vs.z +
                v2.x * vs.y * v1.z +
                vs.x * v1.y * v2.z -
                vs.x * v2.y * v1.z -
                vs.y * v2.z * v1.x -
                vs.z * v2.x * v1.y) * invDenominator;

            if(t < EPS || t > r.max_t) {
                return false;
            }

            float beta = (vs.x * v2.y * v3.z +
                    v2.x * v3.y * vs.z +
                    v3.x * vs.y * v2.z -
                    v3.x * v2.y * vs.z -
                    v3.y * v2.z * vs.x -
                    v3.z * v2.x * vs.y) * invDenominator;


            float gamma = (v1.x * vs.y * v3.z +
                    vs.x * v3.y * v1.z +
                    v3.x * v1.y * vs.z -
                    v3.x * vs.y * v1.z -
                    v3.y * vs.z * v1.x -
                    v3.z * vs.x * v1.y) * invDenominator;

            float alpha = 1.0 - (gamma + beta);
            if(alpha<0 || alpha > 1 || beta<0 || beta>1 || gamma<0 || gamma>1){
                return false;
            }

            Vec2 uv(0.0,0.0);

            Vec3 normal;

            {
                Vec3 n1 = mesh->get_normal(tri[0]);
                Vec3 n2 = mesh->get_normal(tri[1]);
                Vec3 n3 = mesh->get_normal(tri[2]);
                normal = normalized(n1*alpha + n2*beta + n3*gamma);
            }

            if(mesh->has_uv()){
                Vec2 uv1 = mesh->get_uv(tri[0]);
                Vec2 uv2 = mesh->get_uv(tri[1]);
                Vec2 uv3 = mesh->get_uv(tri[2]);
                uv = uv1*alpha + uv2*beta + uv3*gamma;
            }

            Vec3 bitangent(1.0,0.0,0.0);

            if(mesh->use_bitangents()){
                Vec3 t1 = mesh->get_bitangent(tri[0]);
                Vec3 t2 = mesh->get_bitangent(tri[1]);
                Vec3 t3 = mesh->get_bitangent(tri[2]);
                bitangent = normalized(t1*alpha+t2*beta+t3*gamma);
            }

            r.max_t = t;

            if(dot(normal, r.d) >0.0){
                normal = -normal;
            }


            it->normal = normal;
            it->bitangent = bitangent;
            it->tangent = normalized(cross(bitangent,normal));
            it->hitpoint = r.o+r.d*t;
            it->texture_coords = Vec3(uv.x,uv.y,0.0);
            return true;

        }
        bool intersect_any(Ray& r) const {
            Vec3 a,b,c;
            std::array<int,3> tri =  mesh->get_triangle(id);
            a = mesh->get_vertex(tri[0]);
            b = mesh->get_vertex(tri[1]);
            c = mesh->get_vertex(tri[2]);

            //p = a + e1*x+e2*y;
            //p = o +d*t

            //o + d*t = a + e1*x+e2*y
            //o - a = -d*t + e1*x + e2*y

            /*
                       [x]
            [e1 e2 -d] [y] = [o - a]
                       [t]
            */
            Vec3 v1 = b-a;
            Vec3 v2 = c-a;
            Vec3 v3 = -r.d;
            Vec3 vs = (Vec3) (r.o - a);

            float invDenominator = 1.0/(v1.x * v2.y * v3.z +
                                    v2.x * v3.y * v1.z +
                                    v3.x * v1.y * v2.z -
                                    v3.x * v2.y * v1.z -
                                    v3.y * v2.z * v1.x -
                                    v3.z * v2.x * v1.y);

            float t = (v1.x * v2.y * vs.z +
                v2.x * vs.y * v1.z +
                vs.x * v1.y * v2.z -
                vs.x * v2.y * v1.z -
                vs.y * v2.z * v1.x -
                vs.z * v2.x * v1.y) * invDenominator;

            if(t < EPS || t > r.max_t) {
                return false;
            }

            float beta = (vs.x * v2.y * v3.z +
                    v2.x * v3.y * vs.z +
                    v3.x * vs.y * v2.z -
                    v3.x * v2.y * vs.z -
                    v3.y * v2.z * vs.x -
                    v3.z * v2.x * vs.y) * invDenominator;


            float gamma = (v1.x * vs.y * v3.z +
                    vs.x * v3.y * v1.z +
                    v3.x * v1.y * vs.z -
                    v3.x * vs.y * v1.z -
                    v3.y * vs.z * v1.x -
                    v3.z * vs.x * v1.y) * invDenominator;

            float alpha = 1.0 - (gamma + beta);
            if(alpha<0 || alpha > 1 || beta<0 || beta>1 || gamma<0 || gamma>1){
                return false;
            }
            return true;

        }

        AABB get_bounds() const {

            std::array<int,3> tri =  mesh->get_triangle(id);
            Vec3 a,b,c;
            a = mesh->get_vertex(tri[0]);
            b = mesh->get_vertex(tri[1]);
            c = mesh->get_vertex(tri[2]);

            return AABB(
                minF(minF(a,b),c), maxF(maxF(a,b),c)
            );


            //return aabb; let's try to not store the bounds
        }

        float area() const {

            //Todo: store area and aabb?
            Vec3 a,b,c;
            std::array<int,3> tri =  mesh->get_triangle(id);
            a = mesh->get_vertex(tri[0]);
            b = mesh->get_vertex(tri[1]);
            c = mesh->get_vertex(tri[2]);

            return length(cross(b-a,c-a))*0.5;
        }

        Intersection sample(Sampler& s) const{
            Intersection it;

            float r1 = s.sample();
            float r2 = s.sample();

            float alpha = r1;
            float beta = r2;

            if(alpha + beta > 1.0) {
                alpha = 1.0 - alpha;
                beta  = 1.0 - beta;
            }

            float gamma = 1.0 - alpha - beta;
            Vec2 uv(0.0,0.0);


            Vec3 a,b,c;
            std::array<int,3> tri =  mesh->get_triangle(id);
            a = mesh->get_vertex(tri[0]);
            b = mesh->get_vertex(tri[1]);
            c = mesh->get_vertex(tri[2]);

            Vec3 pos = a*alpha + b*beta + c*gamma;


            Vec3 normal;

            {
                Vec3 n1 = mesh->get_normal(tri[0]);
                Vec3 n2 = mesh->get_normal(tri[1]);
                Vec3 n3 = mesh->get_normal(tri[2]);
                normal = normalized(n1*alpha + n2*beta + n3*gamma);
            }

            if(mesh->has_uv()){
                Vec2 uv1 = mesh->get_uv(tri[0]);
                Vec2 uv2 = mesh->get_uv(tri[1]);
                Vec2 uv3 = mesh->get_uv(tri[2]);
                uv = uv1*alpha + uv2*beta + uv2*gamma;
            }

            Vec3 bitangent(1.0,0.0,0.0);

            if(mesh->use_bitangents()){
                Vec3 t1 = mesh->get_bitangent(tri[0]);
                Vec3 t2 = mesh->get_bitangent(tri[1]);
                Vec3 t3 = mesh->get_bitangent(tri[2]);
                bitangent = normalized(t1*alpha+t2*beta+t3*gamma);
            }

            it.normal = normal;
            it.bitangent = bitangent;
            it.tangent = normalized(cross(bitangent,normal));
            it.hitpoint = pos;
            it.texture_coords = Vec3(uv.x,uv.y,0.0);
            return it;
        }

    };
}

#endif