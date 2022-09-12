#ifndef CPPPT_CONSTANT_TEXTURE
#define CPPPT_CONSTANT_TEXTURE

namespace cpppt{



template <typename T>
class ConstantTexture: public Texture {
    T constant;
public:
    ConstantTexture(const T& ct): constant(ct) {}


    Vec3 sample(const Vec3& uv) const {
        return Vec3(constant);
    }
};

}

#endif