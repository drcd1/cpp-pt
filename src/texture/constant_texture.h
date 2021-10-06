#ifndef CPPPT_CONSTANT_TEXTURE
#define CPPPT_CONSTANT_TEXTURE

namespace cpppt{



template <typename T>
class ConstantTexture: public Texture<T> {
    T constant;
public:
    ConstantTexture(const T& ct): constant(ct) {}


    T sample(const Vec3& uv) const {
        return constant;
    }
};

}

#endif