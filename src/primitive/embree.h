#ifndef CPPPT_EMBREE
#define CPPPT_EMBREE

#include <primitive/ray.h>
#include <primitive/primitive.h>
#include <primitive/simple_group.h>
#include <shape/triangle.h>
#include <primitive/aabb.h>
#include <scene.h>

#include <embree4/rtcore.h>
#include <memory>


namespace cpppt{


class Embree: public Primitive {
private:
    RTCDevice device;
    RTCScene scene;
    RTCGeometry geom;
    AABB bounds;
    std::vector<unsigned int> mesh;
    std::vector<unsigned int> nt_cum_sum;
    const SceneData* scene_data;
    std::vector<std::shared_ptr<PrimitiveLeaf>> primitives;

public:
    ~Embree(){
        rtcReleaseScene(scene);
        rtcReleaseDevice(device);
    }

    Embree(const SceneData& sd,
    const std::vector<std::shared_ptr<PrimitiveLeaf>>&& prims)

    :
        device(rtcNewDevice(NULL)),
        scene(rtcNewScene(device)),
        geom(rtcNewGeometry(device,RTC_GEOMETRY_TYPE_TRIANGLE)),
        scene_data(&sd),
        primitives(prims)
    {
        if(!device){
            throw std::runtime_error("Can't init device");
        }

        if(!scene){
            throw std::runtime_error("Can't init scene");
        }

        int tv = 0;
        int tf = 0;

        int actTV = 0;
        int actTF = 0;

        nt_cum_sum.push_back(0);
        for(int i = 0; i<sd.meshes.size(); i++){
            tv+=sd.meshes.at(i)->n_vertices();
            tf+=sd.meshes.at(i)->n_triangles();
            nt_cum_sum.push_back(tf);
        }
        mesh = std::vector<unsigned int>(tf);


        float* vb = (float*) rtcSetNewGeometryBuffer(
            geom,
            RTC_BUFFER_TYPE_VERTEX,
            0,
            RTC_FORMAT_FLOAT3,
            3*sizeof(float),
            tv
        );

        unsigned* ib = (unsigned*) rtcSetNewGeometryBuffer(
            geom,
            RTC_BUFFER_TYPE_INDEX,
            0,
            RTC_FORMAT_UINT3,
            3*sizeof(unsigned),
            tf
        );
        actTV = tv;
        actTF = tf;

        tv = 0;
        tf = 0;


        for(int i = 0; i<sd.meshes.size(); i++){
            for(int j = 0; j<sd.meshes.at(i)->n_vertices();j++){
                const auto& v = sd.meshes.at(i)->get_vertex(j);
                vb[(j+tv)*3]=v.x;
                vb[(j+tv)*3+1]=v.y;
                vb[(j+tv)*3+2]=v.z;
                if((j+tv)*3+2>=actTV*3){
                    std::cout<<"what="<<std::endl;
                }
            }


            for(int j = 0; j<sd.meshes.at(i)->n_triangles();j++){
                const auto& v = sd.meshes.at(i)->get_triangle(j);
                ib[(j+tf)*3]=v[0]+tv;
                ib[(j+tf)*3+1]=v[1]+tv;
                ib[(j+tf)*3+2]=v[2]+tv;
                if((j+tf)*3+2>=actTF*3){
                    std::cout<<"what="<<std::endl;
                }
                mesh.at(j+tf) = i;

            }

            tv+=sd.meshes.at(i)->n_vertices();
            tf+=sd.meshes.at(i)->n_triangles();
            /*
            for(int j = 0; j<sd->meshes.at(i)->n_triangles(); j++){
                auto helper = std::make_shared<PrimitiveLeaf>(
                    std::make_shared<Triangle>(sd->meshes.at(i).get(),j),
                    sd->materials.at(i)
                );

                bvh->add(helper);

                if(sd->materials.at(i)->is_emissive()){
                    ld.lg->add(std::make_shared<ShapeLight>(helper));
                }

            }
            */
        }
    /*
        vb[0] = 0.0f;
        vb[1] = 0.0f;
        vb[2] = 0.0f;


        vb[3] = 0.0f;
        vb[4] = 1.0f;
        vb[5] = 1.0f;


        vb[6] = 0.0f;
        vb[7] = 1.0f;
        vb[8] = 0.0f;

        ib[0] = 0;
        ib[1] = 1;
        ib[2] = 2;
    */

        rtcCommitGeometry(geom);
        rtcAttachGeometry(scene, geom);
        rtcReleaseGeometry(geom);
        rtcCommitScene(scene);


        RTCBounds b;
        rtcGetSceneBounds(scene,&b);
        bounds = AABB(
            {b.lower_x, b.lower_y, b.lower_z},
            {b.upper_x, b.upper_y, b.upper_z}
        );
    }

    bool intersect(Ray& r, Intersection* is) const {
        RTCRayHit rayhit = create_ray_hit(r);

        rtcIntersect1(scene, &rayhit);
        if(rayhit.hit.geomID ==RTC_INVALID_GEOMETRY_ID){
            return false;
        } else {
            unsigned int t_id = rayhit.hit.primID;
            auto tri = Triangle((scene_data->meshes.at(mesh.at(t_id)).get()),
                                t_id-nt_cum_sum.at(mesh.at(t_id)));
            r.max_t = rayhit.ray.tfar;

            is->hitpoint = r.o+r.d*r.max_t;
            Vec2 bary = {rayhit.hit.u, rayhit.hit.v};
            //TODO: make sure it's correct
            tri.fill_intersection_from_barycentric({1.0f-bary.x-bary.y,bary.x,bary.y},is);
            //so far, no material reuse
            //TODO: fixme
            is->primitive = primitives.at(t_id).get();
            is->material = (scene_data->materials.at(mesh.at(t_id)));
            return true;
        }
    }
    bool intersect_any(Ray& r) const {

        RTCRay ray = create_ray(r);

        rtcOccluded1(scene, &ray);

        //if hit, tfar is -inf
        return ray.tfar < -10.0;

    }

    AABB get_bounds() const {
        return bounds;
    }
private:
    RTCRay create_ray(Ray& r) const {
        RTCRay ret;

        ret.org_x = r.o.x;
        ret.org_y = r.o.y;
        ret.org_z = r.o.z;

        ret.dir_x = r.d.x;
        ret.dir_y = r.d.y;
        ret.dir_z = r.d.z;
        ret.tnear  = 0.f;
        ret.tfar   = r.max_t;


        ret.mask = -1;
        ret.flags = 0;

        return ret;

    }
    RTCRayHit create_ray_hit(Ray& r) const {
        RTCRayHit ret;
        ret.ray = create_ray(r);
        ret.hit.geomID = RTC_INVALID_GEOMETRY_ID;
        ret.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
        return ret;

    }
};


}

#endif