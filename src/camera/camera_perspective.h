#ifndef CPPPT_CAMERA_PERSPECTIVE
#define CPPPT_CAMERA_PERSPECTIVE

#include <camera/camera.h>
#include <math/math.h>
#include <iostream>

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


    Ray get_ray(Vec2 xy) const {
        Vec3 direction;
        direction.z = -1.0;
        direction.y = xy.y*tan_fovy;
        direction.x = xy.x*tan_fovy*aspect_ratio;
    
        direction = normalized(coords*direction);

        return Ray(origin, direction);
    };
    RgbImage& get_image()  {
        return image;
    }
};
}
#endif