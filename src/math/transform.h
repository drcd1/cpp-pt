#ifndef CPPPT_TRANSFORM
#define CPPPT_TRANSFORM

#include <math/math.h>

namespace cpppt{
    

class Transform {
    private:
        Mat4 mat;
        Mat4 inv_mat;
    public:
        Transform(const Mat4& mat, const Mat4& inv_mat):mat(mat),inv_mat(inv_mat){}
        Transform compose(const Transform& t){
            return Transform(mat*t.mat, t.inv_mat*inv_mat);
        }
        Transform inverse(const Transform& t){
            return Transform(inv_mat, mat);
        }

        static Transform id(){
            return Transform(Mat4::id(), Mat4::id());
        }
        static Transform translate(Vec3 v){
            return Transform(Mat4::translate(v),Mat4::translate(-v));
        }
        static Transform rotateX(float angle){
            return Transform(Mat4::rotateX(angle),Mat4::rotateX(-angle));
        }
        static Transform rotateY(float angle){
            return Transform(Mat4::rotateY(angle),Mat4::rotateY(-angle));
        }
        static Transform rotateZ(float angle){
            return Transform(Mat4::rotateZ(angle),Mat4::rotateZ(-angle));
        }
};



}

#endif