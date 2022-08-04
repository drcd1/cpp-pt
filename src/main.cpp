

#include <camera/camera_perspective.h>
#include <renderer/dummy_renderer.h>
#include <renderer/Pathtracer.h>
#include <primitive/simple_group.h>
#include <primitive/bvh.h>
#include <shape/mesh.h>
#include <shape/triangle.h>
#include <primitive/primitive_leaf.h>
#include <bxdf/diffuse_bxdf.h>
#include <bxdf/emissive_bxdf.h>
#include <texture/constant_texture.h>
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
    std::cout<<aabb.intersectAny(a)<<std::endl;
    std::cout<<aabb.intersectAny(b)<<std::endl;
    std::cout<<aabb.intersectAny(c)<<std::endl;
    std::cout<<aabb.intersectAny(d)<<std::endl;
    std::cout<<aabb.intersectAny(e)<<std::endl;
    std::cout<<aabb.intersectAny(f)<<std::endl;


}

int main(int argc, char** argv){

    Vec2i resolution{256,256};
    int samples = 2;

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

    std::shared_ptr<DiffuseBxDF> floor_mat = std::make_shared<DiffuseBxDF>( std::make_shared<ConstantTexture<Vec3>>(Vec3(0.8,0.8,0.8)));

    std::shared_ptr<DiffuseBxDF> green_mat = std::make_shared<DiffuseBxDF>( std::make_shared<ConstantTexture<Vec3>>(Vec3(0.05,0.8,0.05)));

    std::shared_ptr<DiffuseBxDF> red_mat = std::make_shared<DiffuseBxDF>( std::make_shared<ConstantTexture<Vec3>>(Vec3(0.8,0.05,0.05)));

    std::shared_ptr<EmissiveBxDF> emit = std::make_shared<EmissiveBxDF>( std::make_shared<ConstantTexture<Vec3>>(Vec3(0.6,0.6,0.6)),6.0f);
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
                                        floor_mat
                                    )                                                    
    );



    std::vector<Vec3> verts;
    std::vector<std::array<int,3>> tris;
    int grid_res = 100;

    for(int i = 0; i<grid_res; i++){
        for(int j = 0; j<grid_res; j++){
            float x = float(i)*(1.0/float(grid_res));
            float y = float(j)*(1.0/float(grid_res));
            verts.push_back(Vec3(x*0.5,(cos(x*2.5)*cos(y*5.0))*0.4,y*0.5));
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
    }

    //g.add(std::make_shared<PrimitiveLeaf>(std::make_shared<Triangle>(&m, 0),green_mat));
    //g.add(std::make_shared<PrimitiveLeaf>(std::make_shared<Triangle>(&m, 1),green_mat));


    g.add(std::make_shared<PrimitiveLeaf>(std::make_shared<Sphere>(Vec3(-0.2,-0.5,-0.8),0.2f),
                                        floor_mat
                                    )                                                    
    );

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
                                        floor_mat
                                    )                                                    
    );

    g.add(std::make_shared<PrimitiveLeaf>(std::make_shared<Sphere>(Vec3(0.0,-1001.0,0.),1000.f),
                                        floor_mat
                                    )                                                    
    );

    g.add(std::make_shared<PrimitiveLeaf>(std::make_shared<Sphere>(Vec3(0.6,0.0,7.05),6.06f),
                                        emit
                                    )                                                    
    );
    
    
    
    g.build();

    Scene scene;
    scene.camera = &camera;
    scene.primitive = &g;
    
    Pathtracer renderer(20);
    //DummyRenderer renderer;


    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    
    auto t1 = high_resolution_clock::now();

    renderer.render(scene, "lol.png");

    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;
    
    std::cout << ms_double.count()*0.001 << "s\n";
    test();

}


/*int main(){

    let floor_mat = &(DiffuseMaterial::new(Vec3::xyz(0.8,0.8,0.8))) as & dyn Material;
    let green_mat =&(DiffuseMaterial{albedo: Vec3::xyz(0.05,0.8,0.05)}) as & dyn Material;
    let red_mat = &(DiffuseMaterial{albedo: Vec3::xyz(0.8,0.05,0.05)}) as & dyn Material;
    let emit = &(EmissionMaterial{light: Vec3::xyz(1.0,1.0,1.0), intensity: 6.0}) as & dyn Material;
    
    let sg =&mut SimpleGroup::new();



    sg.add(Box::new(Sphere{
                        o: Vec3{x:0.3,y:0.0,z:0.1},
                        r: 0.2,
                        mat: floor_mat
                    })
            );

    sg.add(Box::new(Sphere{
                o: Vec3{x:-0.2,y:-0.5,z:-0.8},
                r: 0.2,                
                mat: floor_mat
            })
    );


    sg.add(Box::new(Sphere{
                o: Vec3{x:1001.0,y:0.0,z:0.0},
                r: 1000.,
                
                mat: green_mat
            })
    );

    sg.add(Box::new(Sphere{
            o: Vec3{x:-1001.,y:0.0,z:0.0},
            r: 1000.0,
            mat: red_mat
        })
    );
    sg.add(Box::new(Sphere{
            o: Vec3{x:0.0,y:0.0,z:1001.0},
            r: 1000.0,
                
            mat: floor_mat
        })
    );
    sg.add(Box::new(Sphere{
            o: Vec3{x:0.0,y:0.0,z:-1001.0},
            r: 1000.0,
                
            mat: floor_mat
        })
    );

    sg.add(Box::new(Sphere{
            o: Vec3{x:0.0,y:-1001.0,z:0.0},
            r: 1000.0,
                
            mat: floor_mat
        })
    );

    sg.add(Box::new(Sphere{
        o: Vec3{x:0.6,y:0.0,z:7.05},
        r: 6.06,
        mat: emit
    })
);
    

    let mut c = CameraPerspective::new(
        512,512,
        1.0,
        Vec3{x:0.0,y:1.0,z:0.0},
        Vec3{x:0.0,y:-1.0,z:0.0},
        Vec3{x:0.0,y:0.0,z:1.0},
    );


    let mut s = Scene{
        primitive: sg,
        camera: &mut c,
        lights: Vec::new()
    };

    let renderer = PtRenderer::new(50);
    renderer.render(&mut s, &String::from("lol3.png"));
}*/