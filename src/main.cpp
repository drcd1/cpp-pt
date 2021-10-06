

#include <camera/camera_perspective.h>
#include <renderer/dummy_renderer.h>
#include <renderer/Pathtracer.h>
#include <primitive/simple_group.h>
#include <primitive/primitive_leaf.h>
#include <bxdf/diffuse_bxdf.h>
#include <bxdf/emissive_bxdf.h>
#include <texture/constant_texture.h>
#include <shape/shape.h>
#include <scene.h>


using namespace cpppt;
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


    SimpleGroup g;

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
    Scene scene;
    scene.camera = &camera;
    scene.primitive = &g;
    
    Pathtracer renderer(samples);
    renderer.render(scene, "lol.png");
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