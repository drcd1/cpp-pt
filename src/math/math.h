#ifndef CPPPT_MATH
#define CPPPT_MATH


#define _USE_MATH_DEFINES
#include <cmath>

#include <stdexcept>
#include <iostream>

#define EPS 1e-4


namespace cpppt{


template <typename T>
struct Vector3{
    T x,y,z;

    Vector3(): x(0.0),y(0.0),z(0.0){}
    Vector3(T x, T y, T z): x(x),y(y),z(z){}
    Vector3(T x):x(x),y(x),z(x){}

    Vector3<T> operator+(const Vector3<T>& a) const {
        return {x+a.x,y+a.y,z+a.z};
    }

    Vector3<T> operator-(const Vector3<T>& a) const {
        return {x-a.x,y-a.y,z-a.z};
    }

    Vector3<T> operator/(const Vector3<T>& a) const {
        return {x/a.x,y/a.y,z/a.z};
    }

    Vector3<T> operator*(const Vector3<T>& a) const {
        return {x*a.x,y*a.y,z*a.z};
    }

    Vector3<T> operator/(T a) const {
        return (*this)*(1.0/a);
    }
    Vector3<T> operator*(T a) const {
        return {x*a,y*a,z*a};
    }

};

template <typename T>
T dot(const Vector3<T>& a, const Vector3<T>& b){
    return a.x*b.x+a.y*b.y+a.z*b.z;
}

template <typename T>
T lensqr(const Vector3<T>& a){
    return a.x*a.x+a.y*a.y+a.z*a.z;
}

template <typename T>
T len(const Vector3<T>& a){
    return sqrt(lensqr(a));
}

template <typename T>
Vector3<T> cross(const Vector3<T>& a, const Vector3<T>& b){
    return Vector3<T>(a.y*b.z-a.z*b.y,
                a.z*b.x-a.x*b.z,
                a.x*b.y-a.y*b.x
                );
}

template <typename T>
Vector3<T> normalized(const Vector3<T>& a){
    return a/len(a);
}




template <typename T>
struct Vector2{
    T x,y;

    Vector2(): x(0.0),y(0.0){}
    Vector2(T x, T y): x(x),y(y){}
    Vector2(T x):x(x),y(x){}

    Vector2<T> operator+(const Vector2<T>& a) const {
        return {x+a.x,y+a.y};
    }

    Vector2<T> operator-(const Vector2<T>& a) const {
        return {x-a.x,y-a.y};
    }

    Vector2<T> operator/(const Vector2<T>& a) const {
        return {x/a.x,y/a.y};
    }

    Vector2<T> operator*(const Vector2<T>& a) const {
        return {x*a.x,y*a.y};
    }

    Vector2<T> operator/(T a) const {
        return (*this)*(1.0/a);
    }
    Vector2<T> operator*(T a) const {
        return {x*a,y*a};
    }

};

template <typename T>
T dot(const Vector2<T>& a, const Vector2<T>& b){
    return a.x*b.x+a.y*b.y;
}

template <typename T>
T lensqr(const Vector2<T>& a){
    return a.x*a.x+a.y*a.y;
}

template <typename T>
T len(const Vector2<T>& a){
    return sqrt(lensqr(a));
}


template <typename T>
Vector2<T> normalized(const Vector2<T>& a){
    return a/len(a);
}

typedef Vector3<float> Vec3;
typedef Vector3<int> Vec3i;

typedef Vector2<float> Vec2;
typedef Vector2<int> Vec2i;



struct Mat3{
    float values[9];
    Mat3(){
        values[0] = 0.;
        values[1] = 0.;
        values[2] = 0.;
        values[3] = 0.;
        values[4] = 0.;
        values[5] = 0.;
        values[6] = 0.;
        values[7] = 0.;
        values[8] = 0.;
    }

    Mat3(const Vec3& x,const Vec3& y, const Vec3& z){
        values[0] = x.x;
        values[1] = y.x;
        values[2] = z.x;
        values[3] = x.y;
        values[4] = y.y;
        values[5] = z.y;
        values[6] = x.z;
        values[7] = y.z;
        values[8] = z.z;
    }

    Mat3(float a, float b, float c, float d, float e, float f, float h, float g, float i){
        values[0] = a;
        values[1] = b;
        values[2] = c;
        values[3] = d;
        values[4] = e;
        values[5] = f;
        values[6] = g;
        values[7] = h;
        values[8] = i;
    }
    static Mat3 Mat3::id(){
        return {
                1.0,0.0,0.0,
                0.0,1.0,0.0,
                0.0,0.0,1.0
                };
    }

    static Mat3 Mat3::zeros(){
        return {0.0,0.0,0.0,
                0.0,0.0,0.0,
                0.0,0.0,0.0
                };
    }

    Mat3 operator+(const Mat3& a) const {
        return {
            values[0]+a.values[0],
            values[1]+a.values[1],
            values[2]+a.values[2],
            values[3]+a.values[3],
            values[4]+a.values[4],
            values[5]+a.values[5],
            values[6]+a.values[6],
            values[7]+a.values[7],
            values[8]+a.values[8]        
        };
    }

    Mat3 operator-(const Mat3& a) const {
        return {
            values[0]-a.values[0],
            values[1]-a.values[1],
            values[2]-a.values[2],
            values[3]-a.values[3],
            values[4]-a.values[4],
            values[5]-a.values[5],
            values[6]-a.values[6],
            values[7]-a.values[7],
            values[8]-a.values[8]        
        };
    }
    Mat3 operator*(const Mat3& a) const {
        return {
                values[0]*a.values[0]+values[1]*a.values[3]+values[2]*a.values[6],
                values[0]*a.values[1]+values[1]*a.values[4]+values[2]*a.values[7],
                values[0]*a.values[2]+values[1]*a.values[5]+values[2]*a.values[8],
                
                values[3]*a.values[0]+values[4]*a.values[3]+values[5]*a.values[6],
                values[3]*a.values[1]+values[4]*a.values[4]+values[5]*a.values[7],
                values[3]*a.values[2]+values[4]*a.values[5]+values[5]*a.values[8],
                
                values[6]*a.values[0]+values[7]*a.values[3]+values[8]*a.values[6],
                values[6]*a.values[1]+values[7]*a.values[4]+values[8]*a.values[7],
                values[6]*a.values[2]+values[7]*a.values[5]+values[8]*a.values[8]
        };
    }

    Vec3 operator*(const Vec3& a) const {
        return {
                values[0]*a.x+values[1]*a.y+values[2]*a.z,
                values[3]*a.x+values[4]*a.y+values[5]*a.z,
                values[6]*a.x+values[7]*a.y+values[8]*a.z
        };
    }

    Mat3 operator*(float a) const {
        return {
            values[0]*a,
            values[1]*a,
            values[2]*a,
            values[3]*a,
            values[4]*a,
            values[5]*a,
            values[6]*a,
            values[7]*a,
            values[8]*a
        };
    }

    Mat3 operator/(float a) const {
        return (*this)/a;
    }


    Mat3 inverse() const {
        float a = values[0];
        float b = values[1];
        float c = values[2];
        float d = values[3];   // a b c
        float e = values[4];   // d e f
        float f = values[5];   // g h i
        float g = values[6];
        float h = values[7];
        float i = values[8];

        float det = a*e*i+d*h*c+b*f*g - c*e*g-d*b*i-a*f*h;

        
        return Mat3(
                (e*i-f*h), -(b*i-c*h),  (b*f-c*e),
               -(d*i-f*g),  (a*i-c*g), -(a*f-c*d),
                (d*h-e*g), -(a*h-b*g),  (a*e-b*d)
        )/det;
    }
    Mat3 transpose() const {
        return {
                values[0],
                values[3],
                values[6],

                
                values[1],
                values[4],
                values[7],

                
                values[2],
                values[5],
                values[8],    
        };
    }


};  

bool solve_quadratic(float a, float b, float c, float* t1, float* t2){
    float det = b*b-4.0*a*c;
    if(det<0.0) {
        return false;
    }

    det = sqrt(det);
    float div = 1.0/(2.0*a);
    *t1 = (-b-det)*div;
    *t2 = (-b+det)*div;
    return true;
}

Vec3 sample_hemisphere_cos(float r1, float r2){
    float theta = r1*2.0*M_PI;
    float sin_phi =sqrt(r2);
    float cos_phi = sqrt((1.0-sin_phi*sin_phi));

    return Vec3(cos(theta)*sin_phi,sin(theta)*sin_phi,cos_phi);
}


//gets an orthogonal system -- function inpired PBRT
void orthogonal(Vec3& v1, Vec3* x, Vec3* y, Vec3* z) {
    float abs_x = fabs(v1.x);
    float abs_y = fabs(v1.y);
    float abs_z = fabs(v1.z);

    Vec3 aux;
    if(abs_x<abs_y) {
        if(abs_x < abs_z){
            aux = Vec3(1.0,0.0,0.0);
        } else {
            aux = Vec3(0.0,0.0,1.0);
        }
    } else {
        if(abs_y<abs_z) {
            aux = Vec3(0.0,1.0,0.0);
        } else {
            aux = Vec3(0.0,0.0,1.0);
        }
    }

    Vec3 aux2 = cross(v1,aux);
    *x = normalized(v1);
    *y = normalized(cross(aux2,v1));
    *z = normalized(aux2);
    return;
}

float linear2srgb(float a){
    if(a>0.0031308){
        return (1.0+0.055)*powf(a,1.0/2.4)-0.055;
    } else {
        return 12.92*a;
    }
}

void print(const Mat3& m) {
    std::cout<<m.values[0]<<" "<<m.values[1]<<" "<<m.values[2]<<std::endl;
    std::cout<<m.values[3]<<" "<<m.values[4]<<" "<<m.values[5]<<std::endl;
    std::cout<<m.values[6]<<" "<<m.values[7]<<" "<<m.values[8]<<std::endl;
}


}


#endif