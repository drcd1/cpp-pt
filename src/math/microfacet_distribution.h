#ifndef CPPPT_MICROFACET_DISTRIBUTION
#define CPPPT_MICROFACET_DISTRIBUTION

#include <math.h>
#include <math/sampler.h>

namespace cpppt{
/*GGX aka Trowbridge Reitz taken verbatim from PBRT 3*/
class MicrofacetDistribution{
private:
    float alphaX, alphaY;

    static void sample11(float cosTheta, float U1, float U2,
                                    float *slope_x, float *slope_y) {
    // special case (normal incidence)
        if (cosTheta > .9999) {
            float r = sqrt(U1 / (1 - U1));
            float phi = 6.28318530718 * U2;
            *slope_x = r * cos(phi);
            *slope_y = r * sin(phi);
            return;
        }

        float sinTheta =
            std::sqrt(std::max((float)0, (float)1 - cosTheta * cosTheta));
        float tanTheta = sinTheta / cosTheta;
        float a = 1. / tanTheta;
        float G1 = 2. / (1 + std::sqrt(1.f + 1.f / (a * a)));

        // sample slope_x
        float A = 2 * U1 / G1 - 1;
        float tmp = 1.f / (A * A - 1.f);
        if (tmp > 1e10) tmp = 1e10;
        float B = tanTheta;
        float D = std::sqrt(
            std::max(float(B * B * tmp * tmp - (A * A - B * B) * tmp), float(0)));
        float slope_x_1 = B * tmp - D;
        float slope_x_2 = B * tmp + D;
        *slope_x = (A < 0 || slope_x_2 > 1.f / tanTheta) ? slope_x_1 : slope_x_2;

        // sample slope_y
        float S;
        if (U2 > 0.5f) {
            S = 1.f;
            U2 = 2.f * (U2 - .5f);
        } else {
            S = -1.f;
            U2 = 2.f * (.5f - U2);
        }
        float z =
            (U2 * (U2 * (U2 * 0.27385f - 0.73369f) + 0.46341f)) /
            (U2 * (U2 * (U2 * 0.093073f + 0.309420f) - 1.000000f) + 0.597999f);
        *slope_y = S * z * std::sqrt(1.f + *slope_x * *slope_x);

        if(std::isinf(*slope_y) || std::isnan(*slope_y)){
            *slope_y = 0.0;
        }
    }

public:

    MicrofacetDistribution(float alphaX, float alphaY): alphaX(alphaX), alphaY(alphaY){}
/*
    Vec3 sample(const Vec3& wo,Sampler& samp) const {
        float r1 = samp.sample();
        float r2 = samp.sample();
        //todo: check https://anderslanglands.github.io/posts/ggx/

        const Vec3 v_h =
            normalized(Vec3(alphaX * wo.x, alphaY * wo.y, wo.z));
        // orthonormal basis
        const float lensq = v_h.x * v_h.x + v_h.y * v_h.y;
        const Vec3 T1 = lensq > 0 ? Vec3(-v_h.y, v_h.x, 0.0f) / sqrtf(lensq)
                            : Vec3(1, 0, 0);
        const Vec3 T2 = cross(v_h, T1);
        // parameterization of projected area
        const float r = sqrtf(r1);
        const float phi = 2.0f * M_PI * r2;
        const float t1 = r * cosf(phi);
        float t2 = r * sinf(phi);
        const float s = 0.5f * (1.0f * v_h.z);
        float arg = 1.0f - t1 * t1;
        arg = arg<0.0?0.0:arg;
        t2 = (1.0f - s) * sqrtf(arg) + s * t2;
        // reprojection onto hemisphere
        arg = 1.0f - t1 * t1 - t2 * t2;
        arg = arg<0.0?0.0:arg;
        const Vec3 n_h =
             T1 *t1 +  T2 * t2 + v_h* sqrtf(arg) ;
        // transform back to ellipsoid
        return normalized(Vec3(alphaX * n_h.x, alphaY * n_h.y, std::max(0.0f, n_h.z)));
    }
*/


    Vec3 sample(const Vec3 &wi, Sampler& s) {
        // 1. stretch wi
        Vec3 wiStretched =
            normalized(Vec3(alphaX * wi.x, alphaY * wi.y, wi.z));

        // 2. simulate P22_{wi}(x_slope, y_slope, 1, 1)
        float slope_x, slope_y;
        float U1 = s.sample();
        float U2 = s.sample();
        sample11(cosTheta(wiStretched), U1, U2, &slope_x, &slope_y);

        // 3. rotate
        float tmp = cosPhi(wiStretched) * slope_x - sinPhi(wiStretched) * slope_y;
        slope_y = sinPhi(wiStretched) * slope_x + cosPhi(wiStretched) * slope_y;
        slope_x = tmp;

        // 4. unstretch
        slope_x = alphaX * slope_x;
        slope_y = alphaY * slope_y;
    //std::cout<<"alphaX "<<alphaX<<" and "<<slope_x<<" and "<<alphaY<<" and "<<slope_y<<std::endl;
        // 5. compute normal
        return normalized(Vec3(-slope_x, -slope_y, 1.));
    }

        //Trowbridge–Reitz
    float d(const Vec3& wh) const {
        float tan2theta = (sinTheta(wh)/cosTheta(wh));
        tan2theta = tan2theta*tan2theta;
        if (std::isinf(tan2theta))
            return 0.;
        float cos4theta = cos2Theta(wh) * cos2Theta(wh);
        float e = (cosPhi(wh)*cosPhi(wh) / (alphaX * alphaX) +
                    sinPhi(wh)*sinPhi(wh) / (alphaY * alphaY)) * tan2theta;
        return 1.0 / (M_PI * alphaX * alphaY * cos4theta * (1 + e) * (1 + e));
    }

    //the normal always points up
    float g(const Vec3& wo, const Vec3& wi) const {
        return 1.0/(1.0+lambda(wo) + lambda(wi));
    }

    float lambda(const Vec3& w) const {
        float absTanTheta = std::abs(cosTheta(w)/sinTheta(w));
        if (std::isinf(absTanTheta)) return 0.;
        float alpha = std::sqrt(cosPhi(w)*cosPhi(w) * alphaX * alphaX +
                               sinPhi(w)*sinPhi(w) * alphaY * alphaY);

        float alpha2Tan2Theta = (alpha * absTanTheta) * (alpha * absTanTheta);
        return (-1.0 + sqrt(1.0 + alpha2Tan2Theta)) / 2.0;
    }

    float g1(const Vec3& w) const {
        //    if (Dot(w, wh) * CosTheta(w) < 0.) return 0.;
        return 1.0 / (1.0 + lambda(w));
    }

    float   pdf(const Vec3& wo, const Vec3& wh) const {

        float ad = d(wh);
        float ag1 = g1(wo);
        float adot = fabs(dot(wo, wh));
        float aother_dot = fabs(cosTheta(wo));
        float ret = ad * ag1 * adot / aother_dot;
        if(std::isnan(ret)){
            std::cout<<"Yet another nan"<<std::endl;
        }


        return ret;
    }

};

}

#endif