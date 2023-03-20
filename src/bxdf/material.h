#ifndef CPPPT_MATERIAL
#define CPPPT_MATERIAL

#include <math/math.h>
#include <memory>


namespace cpppt{

class BxDF;

class Material{
    public:
        virtual std::shared_ptr<BxDF> get_bxdf(const Vec3& uv) const = 0;
        virtual bool is_emissive() const {
            return false;
        }
        virtual Vec3 get_normal_mapped(const Vec3& uv, const Vec3& tangent, const Vec3& bitangent, const Vec3& normal) const =0;

};

}
#endif