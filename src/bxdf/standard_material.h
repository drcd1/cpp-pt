#ifndef CPPPT_STANDARD_MATERIAL
#define CPPPT_STANDARD_MATERIAL

//#include "bxdf/disney_bxdf.h"

#include "bxdf/bxdf.h"
#include "bxdf/diffuse_bxdf.h"
#include "bxdf/disney_bxdf.h"
#include "bxdf/mirror_bxdf.h"
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
         alpha(alpha) {
            //std::cout<<"roughness: "<<roughness->sample(Vec3(0.0)).x<<std::endl;
         }

        std::shared_ptr<BxDF> get_bxdf(const Vec3& uv) const {
           if(metal->sample(uv).x < 0.1)

           return std::make_shared<DiffuseBxDF>(albedo->sample(uv));
           /* return std::make_shared<DisneyBxDF>(
            albedo->sample(uv),
            metal->sample(uv).x,
            roughness->sample(uv).x,
            specular->sample(uv).x,
            transparent->sample(uv).x,
            ior->sample(uv).x

            );*/
            else
            return std::make_shared<GlossyBxDF>(roughness->sample(uv).x);
        }

        //TODO:double check normal decoding
        virtual Vec3 get_normal_mapped(const Vec3& uv, const Vec3& tangent, const Vec3& bitangent, const Vec3& normal) const {
            Vec3 new_normal = normal_map->sample(uv);
            new_normal.x = new_normal.x*2.0-1.0;
            new_normal.y = new_normal.y*2.0-1.0;
            //new_normal.z = new_normal.z*2.0-1.0;
            Mat3 m(tangent,bitangent,normal);
            return normalized(m*new_normal);
        }
};

}
#endif