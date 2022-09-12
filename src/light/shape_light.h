#ifndef CPPPT_SHAPE_LIGHT
#define CPPPT_SHAPE_LIGHT

#include <light/light.h>
#include <primitive/primitive_leaf.h>

namespace cpppt{
class ShapeLight: public Light {
private:
    std::shared_ptr<PrimitiveLeaf> primitive;
public:
    ShapeLight(std::shared_ptr<PrimitiveLeaf> primitive) : primitive(primitive){};

    //for NEE
    LightSample connect_eye_path(Sampler& s, const Intersection& from) const {
        LightSample ls;
        Intersection it;
        it = primitive->get_shape()->sample(s);

        ls.position = it.hitpoint;
        ls.intensity = primitive->get_material()->get_bxdf(it.texture_coords)->emit(
            normalized(from.hitpoint - ls.position),
            it)*fabs(dot(it.normal, normalized(from.hitpoint - ls.position)));

        //ls.intensity = fabs(dot(it.normal, normalized(from.hitpoint - ls.position)))*8.0;

        //TODO: change this to come from sample
        //so far, all samples are equiprobable
        ls.pdf = 1.0/primitive->get_shape()->area();
        ls.normal = it.normal;
        ls.ref = this;
        return ls;

        /*TODO: divide by 2pi ??*/
    }

    float pdf(const Intersection& point, Vec3 wo){
       throw "not implemented";
       return 0;
       // return (1.0/(primitive->shape().area()))*primtive->material().pdf(wo);
    }

    LightPathStart sample(Sampler& s) const {

        Intersection it;
        it = primitive->get_shape()->sample(s);

        Vec3 sample_direction;
        auto bxdf = primitive->get_material()->get_bxdf(it.texture_coords);
        float p = bxdf->emit_sample(s, it, &sample_direction);
        LightPathStart lps;
        lps.direction = sample_direction;
        lps.position = it.hitpoint;
        lps.radiance =  bxdf->emit(
            normalized(lps.direction),
            it);

        lps.pdf = (1.0/primitive->get_shape()->area())*p;
        return lps;

    }


    bool is_delta() const {
        return false;
    }
};
}
#endif