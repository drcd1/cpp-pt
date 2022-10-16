#ifndef CPPPT_RGB_IMAGE
#define CPPPT_RGB_IMAGE

#include <math/math.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb/stb_image_write.h>

#define STB_IMAGE_IMPLEMENTATION
#include <stb/stb_image.h>

#include <vector>
#include <string>


namespace cpppt{

namespace pixel_ops{
inline float char2float(unsigned char a){
    return float(a)/255.0;
}

inline void put2char1(float a, unsigned char* pix){
    if(std::isnan(a)){
        std::cout<<"NAN!"<<std::endl;
    }
    if(a<=0.0)
        a = 0.0;
    if(a>=0.99999){
        a=0.99999;
    }

    *pix = unsigned char(unsigned int(a*255.0));
    return;
}
inline void put2char(float a, unsigned char* pix){
    put2char1(a, pix);
    put2char1(a, pix+1);
    put2char1(a, pix+2);
}

inline void put2char(Vec3 a, unsigned char* pix){

    put2char1(a.x, pix);
    put2char1(a.y, pix+1);
    put2char1(a.z, pix+2);
    return;
}


inline void put2float(float a, float* pix){
    *(pix) = a;
    *(pix+1) = a;
    *(pix+2) = a;
    return;
}

inline void put2float(Vec3 a, float* pix){
    *(pix) = a.x;
    *(pix+1) = a.y;
    *(pix+2) = a.z;
    return;
}
template <typename P>
inline unsigned int pixSize(){
    return 1;
}

template <>
inline unsigned int pixSize<Vec3>(){
    return 3;
}

template <int n>
inline void pix_write(Vec3* t, unsigned char* s){
    (*t).x = char2float(s[0]);
    (*t).y = char2float(s[0]);
    (*t).z = char2float(s[0]);
}

template <>
inline void pix_write<3>(Vec3* t, unsigned char* s){
    (*t).x = char2float(s[0]);
    (*t).y = char2float(s[1]);
    (*t).z = char2float(s[2]);
}

template <int n>
inline void pix_write(float* t, unsigned char* s){
    (*t) = char2float(s[0]);
}


template <int n>
inline void pix_write(Vec3* t, float* s){
    (*t).x = (s[0]);
    (*t).y = (s[0]);
    (*t).z = (s[0]);
}

template <>
inline void pix_write<3>(Vec3* t, float* s){
    (*t).x = (s[0]);
    (*t).y = (s[1]);
    (*t).z = (s[2]);
}

template <int n>
inline void pix_write(float* t, float* s){
    (*t) = (s[0]);
}


inline float linear2srgb(float a){
    if(a>0.0031308){
        return (1.0+0.055)*powf(a,1.0/2.4)-0.055;
    } else {
        return 12.92*a;
    }
}


inline float srgb2linear(float a){
    if(a>0.04045){
        return powf((a + 0.055)/(1.055),2.4);
    } else {
        return a/12.92;
    }
}

inline Vec3 linear2srgb(Vec3 v){
    return Vec3(linear2srgb(v.x),
    linear2srgb(v.y),
    linear2srgb(v.z)
    );
}


inline Vec3 srgb2linear(Vec3 v){
    return Vec3(srgb2linear(v.x),
    srgb2linear(v.y),
    srgb2linear(v.z)
    );
}

}
namespace ImageDef{
enum COLOR_SPACE{
    SRGB, LINEAR
};

enum Border {
    CLAMP, REPEAT,MIRROR, MIRROR_Y_REPEAT_X
};
}

template <typename P>
class Image{

public:

    ImageDef::Border border = ImageDef::REPEAT;

    Image(Vec2i res): pixels(res.x*res.y), res(res){}

    Image(std::string filename, ImageDef::COLOR_SPACE is = ImageDef::SRGB){
        int x, y, n;

        if(stbi_is_hdr(filename.c_str())){
            float *data = stbi_loadf(filename.c_str(), &x, &y, &n, 0);
            pixels.resize(x*y);
            res.x = x;
            res.y = y;
            if(n>=3){
                for(int i = 0; i<pixels.size(); i++){
                    pixel_ops::pix_write<3>(&pixels.at(i), &data[i*n]);
                }
            } else {
                for(int i = 0; i<pixels.size(); i++){
                    pixel_ops::pix_write<1>(&pixels.at(i), &data[i*n]);
                }
            }
            free(data);
        } else {
            unsigned char *data = stbi_load(filename.c_str(), &x, &y, &n, 0);
            pixels.resize(x*y);
            res.x = x;
            res.y = y;
            if(n>=3){
                for(int i = 0; i<pixels.size(); i++){
                    pixel_ops::pix_write<3>(&pixels.at(i), &data[i*n]);
                }
            } else {
                for(int i = 0; i<pixels.size(); i++){
                    pixel_ops::pix_write<1>(&pixels.at(i), &data[i*n]);
                }
            }
            free(data);
        }

        if(is==ImageDef::SRGB){
            for(int i = 0; i<pixels.size(); i++){
                pixels.at(i) = pixel_ops::srgb2linear(pixels.at(i));
            }
        }
    }

    void put_pixel(int i, int j, P pixel){
        pixels.at(i+j*res.x) = pixel;
    }
    P* get_pixel(int i, int j){
        handle_border(&i,&j);
        return &(  pixels.data()[i+j*res.x]);
    }
    void handle_border(int* i, int* j){

        if(*i >= res.x){
            switch(border){
                case ImageDef::CLAMP:
                    *i = res.x;
                break;
                case ImageDef::REPEAT:
                case ImageDef::MIRROR_Y_REPEAT_X:
                    *i = (*i)%res.x;
                break;
                case ImageDef::MIRROR:
                default:
                    int help = *i/res.x;
                    if(help%2 == 0){
                        *i = (*i)%res.x;
                    } else {
                        *i = (-(*i)-1)%res.x + res.x;
                    }
                break;

            }
        } else if (*i<0){
            switch(border){
                case ImageDef::CLAMP:
                    *i = 0;
                break;
                case ImageDef::REPEAT:
                case ImageDef::MIRROR_Y_REPEAT_X:
                    *i = (*i)%res.x + res.x;
                break;
                case ImageDef::MIRROR:
                default:
                    int help = *i/res.x;
                    if(help%2 == 0){
                        *i = (*i)%res.x + res.x;
                    } else {
                        *i = (-(*i)-1)%res.x;
                    }
                break;
            }
        }

        if(*j >= res.y){
            switch(border){
                case ImageDef::CLAMP:
                    *j = res.y;
                break;
                case ImageDef::REPEAT:
                    *j = (*j)%res.y;
                break;
                case ImageDef::MIRROR:
                case ImageDef::MIRROR_Y_REPEAT_X:
                default:
                    int help = *j/res.y;
                    if(help%2 == 0){
                        *j = (*j)%res.y;
                    } else {
                        *j = (-(*j)-1)%res.y + res.y;
                    }
                break;

            }
        } else if (*j<0){
            switch(border){
                case ImageDef::CLAMP:
                    *j = 0;
                break;
                case ImageDef::REPEAT:
                    *j = (*j)%res.y + res.y;
                break;
                case ImageDef::MIRROR:
                case ImageDef::MIRROR_Y_REPEAT_X:
                default:
                    int help = *j/res.y;
                    if(help%2 == 0){
                        *j = (*j)%res.y + res.y;
                    } else {
                        *j = (-(*j)-1)%res.y;
                    }
                break;
            }
        }
    }

    void save(std::string filename, ImageDef::COLOR_SPACE is = ImageDef::SRGB){
        {
            std::vector<unsigned char> data(res.x*res.y*3);

            for(int i = 0; i<pixels.size(); i++){
                auto p = pixel_ops::linear2srgb(pixels.at(i));
                pixel_ops::put2char(p,&data.at(i*3));
            }
            stbi_write_png((filename+".png").c_str(), res.x,res.y,3, (void*)data.data(),res.x*3);

        }
        {
            std::vector<float> data(res.x*res.y*3);
            for(int i = 0; i<pixels.size(); i++){
                pixel_ops::put2float(pixels.at(i),&data.at(i*3));
            }
            stbi_write_hdr((filename+".hdr").c_str(), res.x,res.y,3,data.data());
        }
    }

    Vec2i res;
private:
    std::vector<P> pixels;
};

typedef Image<Vec3> RgbImage;

}
#endif