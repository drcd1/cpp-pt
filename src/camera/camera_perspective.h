#ifndef CPPPT_CAMERA_PERSPECTIVE
#define CPPPT_CAMERA_PERSPECTIVE

#include <camera/camera.h>
#include <math/math.h>
#include <iostream>
#include <algorithm>


namespace cpppt{
class CameraPerspective: public Camera {
private:
    RgbImage image;
    float tan_fovy;
    Vec3 origin;
    Mat3 coords;
    float aspect_ratio;
public:
    CameraPerspective(Vec2i resolution, float tan_fovy, Vec3 origin, Vec3 forward, Vec3 up): image(resolution), tan_fovy(tan_fovy), origin(origin){
        Vec3 x,y,z;
        z = normalized(forward*(-1.0));
        x = normalized(cross(up,z));
        y = normalized(cross(z,x));
        coords = Mat3(x,y,z);
        print(coords);
        aspect_ratio = float(resolution.x)/float(resolution.y);
    }

    float getFovy() const {
        return 2.0*atan(tan_fovy);
    }
    const Vec3& getOrigin() const {
        return origin;
    }

    const Mat3& getCoords() const {
        return coords;
    }

    Ray get_ray(Vec2 xy) const {
        Vec3 direction;
        direction.z = -1.0;
        direction.y = xy.y*tan_fovy;
        direction.x = xy.x*tan_fovy*aspect_ratio;

        direction = normalized(coords*direction);

        return Ray(origin, direction);
    };


    Ray get_ray(Vec2 xy, float& angular_pdf) const {
        Vec3 direction;
        direction.z = -1.0;
        direction.y = xy.y*tan_fovy;
        direction.x = xy.x*tan_fovy*aspect_ratio;

        Vec3 d = direction;
        float r = length(direction);

        direction = normalized(coords*direction);

        float area_pdf = 1.0/(tan_fovy*tan_fovy*aspect_ratio);


        angular_pdf = area_pdf*r*r*r/fabs(d.z);

        return Ray(origin, direction);
    };
    RgbImage& get_image()  {
        return image;
    }



    virtual float pdf(const Vec3& pt, const Vec3& normal) const {
        float p = 1.0/(image.res.x*image.res.y);
        Vec3 dir = coords.transpose()*(pt-origin);
        //TODO: optimize
        float r = length(dir);
        Vec3 n_dir = dir/r;
        p *= dot(dir,coords*(Vec3(0.0,0.0,-1.0)));
        p *= fabs(dot(normal,n_dir))/(r*r);
        return p;

    }


    CameraConnection connect_light_path(Sampler& s, const Intersection& it) const {
        CameraConnection cc;
        cc.pos = origin;

        Vec3 dir = coords.transpose()*(it.hitpoint-cc.pos);
        float n_d_l = length(dir);
        if(n_d_l<EPS || dir.z>=0){
            cc.i=-1;
            cc.j=-1;
            return cc;
        }
        Vec3 n_dir = normalized(dir/n_d_l);
        dir = dir/dir.z * -1.0;
        dir.x =dir.x/(tan_fovy*aspect_ratio);
        dir.y = dir.y/tan_fovy;
        cc.i = int((dir.x*0.5 + 0.5)*image.res.x);
        cc.j = int((dir.y*0.5 + 0.5)*image.res.y);
        if((dir.x*0.5 + 0.5)*image.res.x<0.0 || (dir.y*0.5 + 0.5)*image.res.y<0.0){
            cc.i = -1;
            cc.j = -1;
            return cc;
        }

        cc.j =  image.res.y - cc.j -1;
        if(cc.i>=image.res.x || cc.j >= image.res.y || cc.j<0 || cc.i<0){
            cc.i = -1;
            cc.j = -1;
        }
        //todo: account for x
        //TODO: double check cosine
        Vec3 forwards = coords*Vec3(0.0,0.0,1.0);

        //CAMERA PLANE GOES -1 TO 1, SO AREA FACTOR IS 1/4
        cc.factor = /*std::fabs(dot(it.normal,forwards))/*/1.0/(4.0*std::fabs(n_dir.z*n_dir.z*n_dir.z) *n_d_l*n_d_l *tan_fovy* tan_fovy);

        return cc;


    }
};
}
#endif