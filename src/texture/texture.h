#ifndef CPPPT_TEXTURE
#define CPPPT_TEXTURE

namespace cpppt{



template <typename T>
class Texture {
public:
    enum Filtering {
        POINT, BILINEAR, BICUBIC
    };

    enum Border {
        CLIP, REPEAT,MIRROR, MIRROR_Y_REPEAT_X
    };

    virtual T sample(const Vec3& uv) const = 0;

    static Vec2 handle_border(Border b, const Vec2& uv){
        switch(b){
            case CLIP:
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

};


}

#endif