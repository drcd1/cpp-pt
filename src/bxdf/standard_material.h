#ifndef CPPPT_STANDARD_MATERIAL
#define CPPPT_STANDARD_MATERIAL

//#include "bxdf/disney_bxdf.h"

#include "bxdf/bxdf.h"
#include "bxdf/diffuse_bxdf.h"
#include "texture/texture.h"
#include "texture/constant_texture.h"
#include "math/math.h"
#include "bxdf/material.h"
namespace cpppt{
class StandardMaterial : public Material{
    private:
        std::shared_ptr<Texture> albedo;
        std::shared_ptr<Texture> normal_map;
        std::shared_ptr<Texture> specular;
        std::shared_ptr<Texture> roughness;
        std::shared_ptr<Texture> metal;
        std::shared_ptr<Texture> ior;
        std::shared_ptr<Texture> transparent;
        std::shared_ptr<Texture> alpha;

    public:
        StandardMaterial(
        std::shared_ptr<Texture> albedo,
        std::shared_ptr<Texture> normal_map,
        std::shared_ptr<Texture> specular,
        std::shared_ptr<Texture> roughness,
        std::shared_ptr<Texture> metal,
        std::shared_ptr<Texture> ior,
        std::shared_ptr<Texture> transparent,
        std::shared_ptr<Texture> alpha): albedo(albedo),
        normal_map(normal_map),
         specular(specular),
        roughness(roughness),
         metal(metal),
         ior(ior),
         transparent(transparent),
        alpha(alpha) {}

        std::shared_ptr<BxDF> get_bxdf(Vec3 uv) const {
            return std::make_shared<DiffuseBxDF>(albedo->sample(uv));
        }

};

}
#endif