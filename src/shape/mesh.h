#ifndef CPPPT_MESH
#define CPPPT_MESH

#include <math/math.h>
#include <vector>
#include <array>
namespace cpppt{
class Mesh{
    private:
        std::vector<Vec3> vertices;
        std::vector<Vec2> uvs;
        std::vector<Vec3> normals;

        //bitangents are aligned with V
        //we store bitangents and not tangents because this is what blender does
        //when baking normal maps.
        std::vector<Vec3> bitangents;



        std::vector<std::array<int,3>> triangles;

        void compute_vertex_normals(){
            //todo: check init
            for(auto t: triangles){
                Vec3 a,b,c;
                a=vertices.at(t.at(0));
                b=vertices.at(t.at(1));
                c=vertices.at(t.at(2));

                Vec3 normal = normalized(cross((b-a),(c-a)));
                normals.at(t.at(0)) = normals.at(t.at(0))+normal;
                normals.at(t.at(1)) = normals.at(t.at(1))+normal;
                normals.at(t.at(2)) = normals.at(t.at(2))+normal;
            }

            for(auto& n: normals){
                n = normalized(n);
            }
        }

        void compute_bitangents(){
            if(uvs.size() < 1){
                return;
            }

            for(auto t: triangles){
                Vec3 a = vertices.at(t.at(1))-vertices.at(t.at(0));
                Vec3 b = vertices.at(t.at(2))-vertices.at(t.at(0));

                Vec2 au = uvs.at(t.at(1))-uvs.at(t.at(0));
                Vec2 bu = uvs.at(t.at(2))-uvs.at(t.at(0));

                //if aligned with x
                //Vec3 tangent = normalized(a*bu.y-b*au.y);

                //if aligned with y
                Vec3 tangent = normalized(-a*bu.x + b*au.x)/(au.x*bu.y - au.y*bu.x);
                bitangents.at(t.at(0)) =bitangents.at(t.at(0))+tangent;
                bitangents.at(t.at(1)) =bitangents.at(t.at(1))+tangent;
                bitangents.at(t.at(2)) =bitangents.at(t.at(2))+tangent;
            }

            for(auto& t: bitangents){
                t = normalized(t);
            }



        }


        /*
        static std::vector<Mesh> read_sob(std::string filename){
            //nv nf
            //v
            //10 12 13
            //f
            //tg
            //uv
        }
        */


    public:

        Mesh(
            std::vector<std::array<int,3>>&& triangles,
            std::vector<Vec3>&& vertices,
            std::vector<Vec2>&& uvs = std::vector<Vec2>()
        ):triangles(triangles),vertices(vertices),uvs(uvs)
        {
            normals.resize(vertices.size());
            compute_vertex_normals();
            bitangents.resize(vertices.size());
            compute_bitangents();
        }

        Mesh(int v, int f): triangles(f),vertices(v),normals(v),uvs(v),bitangents(v){}

        void set_vertex_data(int id, Vec3 v, Vec3 n, Vec2 uv, Vec3 bitangent){
            vertices.at(id) = v;
            normals.at(id) = n;
            uvs.at(id) = uv;
            bitangents.at(id) = bitangent;
        }

        void set_face_data(int id, int a, int b, int c){
            triangles.at(id) = {a,b,c};
        }

        Vec3 get_vertex(int id) const {
            return vertices.at(id);
        }
        std::array<int,3> get_triangle(int id) const {
            return triangles.at(id);
        }
        Vec2 get_uv(int id) const {
            return uvs.at(id);
        }
        Vec3 get_bitangent(int id) const {
            return bitangents.at(id);
        }
        Vec3 get_normal(int id) const {
            return normals.at(id);
        }


        bool has_uv() const {
            return uvs.size()>0;
        }

        bool use_bitangents() const {
            return bitangents.size()>0;
        }

        int n_triangles() const {
            return triangles.size();
        }

        int n_vertices() const {
            return vertices.size();
        }

        const auto& getVertices() const {
            return vertices;
        }

        const auto& getUvs() const {
            return uvs;
        }

        const auto& getNormals() const {
            return normals;
        }

        const auto& getBitangents() const {
            return bitangents;
        }
        const auto& getTriangles() const {
            return triangles;
        }



};
}

#endif