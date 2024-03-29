#ifndef CPPPT_EMISSIVE_MATERIAL
#define CPPPT_EMISSIVE_MATERIAL

#include "bxdf/emissive_bxdf.h"

#include "bxdf/bxdf.h"
#include "bxdf/emissive_bxdf.h"
#include "texture/texture.h"
#include "texture/constant_texture.h"
#include "math/math.h"
#include "bxdf/material.h"
namespace cpppt{
class EmissiveMaterial : public Material{
    private:
        std::shared_ptr<Texture> albedo;
        float intensity;
        bool double_sided;

    public:
        EmissiveMaterial(
        std::shared_ptr<Texture> albedo,
        float intensity,bool double_sided = true):
        albedo(albedo),
        intensity(intensity),
        double_sided(double_sided)
        {}

        std::shared_ptr<BxDF> get_bxdf(const Vec3& uv) const {
            return std::make_shared<EmissiveBxDF>(albedo->sample(uv)*intensity, double_sided);
        }

        Vec3 get_normal_mapped(const Vec3& uv, const Vec3& tangent, const Vec3& bitangent, const Vec3& normal) const {
            return normal;
        }

        bool is_emissive() const {
            return true;
        }

};

}
#endif