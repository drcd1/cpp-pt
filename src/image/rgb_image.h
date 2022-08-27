#ifndef CPPPT_RGB_IMAGE
#define CPPPT_RGB_IMAGE

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <math/math.h>

#include <stb/stb_image_write.h>
#include <vector>
#include <string>


namespace cpppt{
char float2char(float a){
    if(a<=0.0)
        a = 0.0;
    if(a>=0.99999){
        a=0.99999;
    }

    return char(unsigned int(a*255.0));
}

class RgbImage{
public:
    RgbImage(Vec2i res): pixels(res.x*res.y), res(res){}
    void put_pixel(int i, int j, Vec3 pixel){
        pixels.at(i+j*res.x) = pixel;
    }
    Vec3* get_pixel(int i, int j){
        return &(  pixels.data()[i+j*res.x]);
    }
    void save(std::string filename){

        std::vector<char> data(res.x*res.y*3);

        for(int i = 0; i<pixels.size(); i++){
            data.at(i*3) = float2char(pixels.at(i).x);
            data.at(i*3+1) = float2char(pixels.at(i).y);
            data.at(i*3+2) = float2char(pixels.at(i).z);
        }

        stbi_write_png(filename.c_str(), res.x,res.y,3, (void*)data.data(),res.x*3);
    }

    Vec2i res;
private:
    std::vector<Vec3> pixels;
};
}
#endif