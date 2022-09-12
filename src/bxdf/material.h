#ifndef CPPPT_MATERIAL
#define CPPPT_MATERIAL

#include "bxdf/bxdf.h"
#include "bxdf/material.h"
namespace cpppt{
class Material{
    public:
        virtual std::shared_ptr<BxDF> get_bxdf(Vec3 uv) const = 0;
        virtual bool is_emissive() const {
            return false;
        }

};

}
#endif