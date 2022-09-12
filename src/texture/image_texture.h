#ifndef CPPPT_IMAGE_TEXTURE
#define CPPPT_IMAGE_TEXTURE

#include <image/rgb_image.h>
#include <sstream>
namespace cpppt{



template <typename T>
class ImageTexture: public Texture {
private:
    std::shared_ptr<Image<T>> im;
public:
    enum Filtering {
        POINT, BILINEAR
    };

    Filtering filtering;

    ImageTexture(std::shared_ptr<Image<T>> im): im(im){
    }

    virtual Vec3 sample(const Vec3& uv) const {
        float x = uv.x*im->res.x;
        float y = (1.0-uv.y)*im->res.y;
        if(filtering == Filtering::POINT){
            return *(im->get_pixel(std::floor(x), std::floor(y)));
        } else /*filtering is bilinear*/{
            x = x-0.5;
            y = y-0.5;
            int i = std::floor(x);
            int j = std::floor(y);
            x = x-float(i);
            y = y-float(j);
            return ((*(im->get_pixel(i,j)))*(1.0-x) + (*(im->get_pixel(i+1,j)))*x )*(1.0-y)+
                    ((*(im->get_pixel(i,j+1)))*(1.0-x) + (*(im->get_pixel(i+1,j+1)))*x)*y;

        }
    }
    /*
    static Vec2 handle_border(Border b, const Vec2& uv){
        switch(b){
            case CLAMP:
                return {clamp(uv.x,0.0,1.0),clamp(uv.y,0.0,1.0)};
            case REPEAT:
                return {fract(uv.x), fract(uv.y)}
            case MIRROR:
                Vec2 uv = {fract(uv.x*0.5)*2.0, fract(uv.y*0.5)*2.0};
                return {uv.x<1.0? uv.x: 2.0 - uv.x, uv.y<1.0? uv.y: 2.0 - uv.y}
            case MIRROR_Y_REPEAT_X:
                Vec2 uv = {fract(uv.x), fract(uv.y*0.5)*2.0};
                return {uv.x, uv.y<1.0? uv.y: 2.0 - uv.y}
            break:
                throw std::runtime_error("Unknown type of border!");
        }
    }
    */

};


}

#endif