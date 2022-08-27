

#include <camera/camera_perspective.h>
#include <renderer/dummy_renderer.h>
#include <renderer/Pathtracer.h>
#include <renderer/Lighttracer.h>
#include <primitive/simple_group.h>
#include <primitive/bvh.h>
#include <shape/mesh.h>
#include <shape/triangle.h>
#include <primitive/primitive_leaf.h>
#include <bxdf/diffuse_bxdf.h>
#include <bxdf/emissive_bxdf.h>
#include <texture/constant_texture.h>
#include <bxdf/mirror_bxdf.h>
#include <bxdf/refraction_bxdf.h>

#include <light/light.h>
#include <light/point_light.h>
#include <light/shape_light.h>
#include <light/light_group.h>

#include <shape/shape.h>
#include <scene.h>
#include <chrono>


using namespace cpppt;

void test(){
    AABB aabb(Vec3(0.0,0.0,0.0),Vec3(1.0,2.0,3.0));
    Ray a({-1.0,-1.0,-1.0},{1.0,1.0,1.0});
    Ray b({0.5,0.5,0.5},{1.0,1.0,1.0});
    Ray c({0.5,0.5,0.5},{-1.0,-1.0,-1.0});


    Ray d({-200.0,-10.0,-0.5},{200.0,10.0,0.5});
    Ray e({-200.0,-10.0,-0.5},{200.0,10.0,1.0});
    Ray f({-200.0,-10.0,-0.5},{200.0,10.0,0.1});
    std::cout<<aabb.intersect_any(a)<<std::endl;
    std::cout<<aabb.intersect_any(b)<<std::endl;
    std::cout<<aabb.intersect_any(c)<<std::endl;
    std::cout<<aabb.intersect_any(d)<<std::endl;
    std::cout<<aabb.intersect_any(e)<<std::endl;
    std::cout<<aabb.intersect_any(f)<<std::endl;


}


int main(int argc, char** argv){

    Vec2i resolution{256,256};
    int samples = 50;

    if(argc >= 3)
        resolution = {atoi(argv[1]), atoi(argv[2])};

    if(argc >= 4)
        samples = atoi(argv[3]);

    Vec3 up = Vec3(0.0,0.0,1.0);
    Vec3 forward = Vec3(0.0,-1.0,0.0);
    Vec3 origin = Vec3(0.0,1.0,0.0);
    CameraPerspective camera(
            resolution,
            1.0,
            origin,
            forward,
            up);


    BVH g;


    Mesh cornell_box(
    {
        {0,1,2},{3,4,5}, //left wall
        {6,7,8},{9,10,11}, //right wall

        {12,13,14},{15,16,17}, //floor
        {18,19,20},{21,22,23}, //ceiling

        {24,25,26}, {27,28,29} //back_wall
    }
    ,
    {


       { 1.0,-1.0,-1.0},
       { 1.0,-1.0, 1.0},
       { 1.0, 1.0,-1.0},


       { 1.0, 1.0,-1.0},
       { 1.0, 1.0, 1.0},
       { 1.0,-1.0, 1.0},

       {-1.0,-1.0,-1.0},
       {-1.0,-1.0, 1.0},
       {-1.0, 1.0,-1.0},

       {-1.0, 1.0,-1.0},
       {-1.0, 1.0, 1.0},
       {-1.0,-1.0, 1.0},


       {-1.0,-1.0,-1.0},
       {-1.0, 1.0,-1.0},
       { 1.0,-1.0,-1.0},

       { 1.0,-1.0,-1.0},
       {-1.0, 1.0,-1.0},
       { 1.0, 1.0,-1.0},


       {-1.0,-1.0,1.0},
       {-1.0, 1.0,1.0},
       { 1.0,-1.0,1.0},

       { 1.0,-1.0,1.0},
       { 1.0, 1.0,1.0},
       {-1.0, 1.0,1.0},

       {-1.0,-1.0,-1.0},
       {-1.0,-1.0, 1.0},
       { 1.0,-1.0,-1.0},

       { 1.0,-1.0,-1.0},
       { 1.0,-1.0, 1.0},
       {-1.0,-1.0, 1.0}

    }
    );

    std::shared_ptr<DiffuseBxDF> floor_mat = std::make_shared<DiffuseBxDF>( std::make_shared<ConstantTexture<Vec3>>(Vec3(0.8,0.8,0.8)));

    std::shared_ptr<DiffuseBxDF> green_mat = std::make_shared<DiffuseBxDF>( std::make_shared<ConstantTexture<Vec3>>(Vec3(0.05,0.8,0.05)));

    std::shared_ptr<DiffuseBxDF> red_mat = std::make_shared<DiffuseBxDF>( std::make_shared<ConstantTexture<Vec3>>(Vec3(0.8,0.05,0.05)));

   // std::shared_ptr<MirrorBxDF> mirror_mat = std::make_shared<MirrorBxDF>();
    std::shared_ptr<MirrorBxDF> mirror_mat = std::make_shared<MirrorBxDF>();
  std::shared_ptr<RefractionBxDF> glass_mat = std::make_shared<RefractionBxDF>(1.4);

/*
    g.add(std::make_shared<PrimitiveLeaf>(std::make_shared<Sphere>(Vec3(-0.2,0.5,0.00),0.1f),
                                            std::make_shared<DiffuseBxDF>(
                                                std::make_shared<ConstantTexture<Vec3>> (Vec3(1.0,0.0,1.0))

                                            )
                                        )
        );

    g.add(std::make_shared<PrimitiveLeaf>(std::make_shared<Sphere>(Vec3(0.00,0.5,0.2),0.1f),
                                     std::make_shared<DiffuseBxDF>(
                                                std::make_shared<ConstantTexture<Vec3>> (Vec3(0.4,0.6,0.1))

                                            )));

    g.add(std::make_shared<PrimitiveLeaf>(std::make_shared<Sphere>(Vec3(0.0,0.5,0.0),0.1f),
                                                    std::make_shared<DiffuseBxDF>(
                                                std::make_shared<ConstantTexture<Vec3>> (Vec3(0.7,0.2,0.1))

                                            )));
*/

    g.add(std::make_shared<PrimitiveLeaf>(std::make_shared<Sphere>(Vec3(0.3,0.0,0.1),0.2f),
                                        mirror_mat)

    );


/*
    std::vector<Vec3> verts;
    std::vector<std::array<int,3>> tris;
    int grid_res = 100;

    for(int i = 0; i<grid_res; i++){
        for(int j = 0; j<grid_res; j++){
            float x = float(i)*(1.0/float(grid_res));
            float y = float(j)*(1.0/float(grid_res));
            verts.push_back(Vec3(x*0.5,(cos(x*5.5)*cos(y*15.0))*0.1,y*0.5));
            if(j<grid_res-1 && i<grid_res-1){
                tris.push_back({i*grid_res+j,(i+1)*grid_res + j,i*grid_res + j+1});
                tris.push_back({i*grid_res+j+1, (i+1)*grid_res + j, (i+1)*grid_res +j+1});
            }
        }
    }




    //Mesh m({ {0,1,2}, {2,1,3}}, {{0.3,0.0,0.1}, {0.3+0.3,0.0,0.1},{0.3,0.2,0.1+0.3}, {0.3+0.3,0.2,0.1+0.3} });
    Mesh m(std::move(tris),std::move(verts));
    for(int i = 0; i<m.n_triangles(); i++){
        //std::cout<<m.get_triangle(i)[0]<<", "<<m.get_triangle(i)[1]<<", "<<m.get_triangle(i)[2]<<std::endl;
        g.add(std::make_shared<PrimitiveLeaf>(std::make_shared<Triangle>(&m, i),green_mat));
    }*/

    //g.add(std::make_shared<PrimitiveLeaf>(std::make_shared<Triangle>(&m, 0),green_mat));
    //g.add(std::make_shared<PrimitiveLeaf>(std::make_shared<Triangle>(&m, 1),green_mat));



    g.add(std::make_shared<PrimitiveLeaf>(std::make_shared<Sphere>(Vec3(-0.2,-0.5,-0.8),0.2f),
                                        glass_mat
                                    )
    );

/*
    g.add(std::make_shared<PrimitiveLeaf>(std::make_shared<Sphere>(Vec3(1001.,0.0,0.0),1000.f),
                                        green_mat
                                    )
    );

    g.add(std::make_shared<PrimitiveLeaf>(std::make_shared<Sphere>(Vec3(-1001.,0.0,0.0),1000.f),
                                        red_mat
                                    )
    );

    g.add(std::make_shared<PrimitiveLeaf>(std::make_shared<Sphere>(Vec3(0.0,0.0,1001.),1000.f),
                                        floor_mat
                                    )
    );

    g.add(std::make_shared<PrimitiveLeaf>(std::make_shared<Sphere>(Vec3(0.0,0.0,-1001.),1000.f),
                                        mirror_mat
                                    )
    );

    g.add(std::make_shared<PrimitiveLeaf>(std::make_shared<Sphere>(Vec3(0.0,-1001.0,0.),1000.f),
                                        floor_mat
                                    )
    );
*/

    g.add(std::make_shared<PrimitiveLeaf>(std::make_shared<Triangle>(&cornell_box,0),
                                        green_mat
                                    )
    );

    g.add(std::make_shared<PrimitiveLeaf>(std::make_shared<Triangle>(&cornell_box,1),
                                        green_mat
                                    )
    );

    g.add(std::make_shared<PrimitiveLeaf>(std::make_shared<Triangle>(&cornell_box,2),
                                        red_mat
                                    )
    );

      g.add(std::make_shared<PrimitiveLeaf>(std::make_shared<Triangle>(&cornell_box,3),
                                        red_mat
                                    )
    );

    g.add(std::make_shared<PrimitiveLeaf>(std::make_shared<Triangle>(&cornell_box,4),
                                        floor_mat
                                    )
    );

    g.add(std::make_shared<PrimitiveLeaf>(std::make_shared<Triangle>(&cornell_box,5),
                                        floor_mat
                                    )
    );

    g.add(std::make_shared<PrimitiveLeaf>(std::make_shared<Triangle>(&cornell_box,6),
                                        floor_mat
                                    )
    );

    g.add(std::make_shared<PrimitiveLeaf>(std::make_shared<Triangle>(&cornell_box,7),
                                        floor_mat
                                    )
    );

    g.add(std::make_shared<PrimitiveLeaf>(std::make_shared<Triangle>(&cornell_box,8),
                                        floor_mat
                                    )
    );

    g.add(std::make_shared<PrimitiveLeaf>(std::make_shared<Triangle>(&cornell_box,9),
                                        floor_mat
                                    )
    );


    //auto light_source = std::make_shared<PrimitiveLeaf>(std::make_shared<Sphere>(Vec3(-0.7,-0.7,-0.45),0.1f),
    //                                   emit);



    std::shared_ptr<EmissiveBxDF> emit =
    std::make_shared<EmissiveBxDF>(
        std::make_shared<ConstantTexture<Vec3>>(
            Vec3(0.6,0.6,0.6)),9.0f);

    LightGroup l;
    //#define AREA_LIGHT
    #ifdef AREA_LIGHT
    auto light_source = std::make_shared<PrimitiveLeaf>(std::make_shared<Sphere>(Vec3(-0.2,-0.2,-0.05),0.1f),
                                      emit);

    g.add(light_source);
    l.add(std::make_shared<ShapeLight>(light_source));
    #else
    l.add(std::make_shared<PointLight>(Vec3(-0.7,-0.4,-0.35), Vec3(0.6,0.6,0.6),9.0f*(4.0*M_PI*0.1*0.1*M_PI)));
    #endif




    g.build();
    l.build();

    Scene scene;
    scene.camera = &camera;
    scene.primitive = &g;
    scene.light = &l;

    Lighttracer renderer(40);
    //DummyRenderer renderer;


    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    auto t1 = high_resolution_clock::now();

    renderer.render(scene, "lighttrace_direct.png");

    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;

    std::cout << ms_double.count()*0.001 << "s\n";
    test();

}
