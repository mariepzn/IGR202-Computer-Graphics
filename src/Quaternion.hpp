
#ifndef _QUATERNION_HPP_
#define _QUATERNION_HPP_
#include "Vector3.hpp"
#include "Matrix3x3.hpp"
#include "math.h"
class Quaternion{
    public :
    Real s;
    Vec3f w;
    Quaternion(Real s, Vec3f w){
        this->s=s;
        this->w=w;
    }
    Quaternion operator/=(Quaternion &d)const{
        return Quaternion(*this) /= d;
    }
    Quaternion operator/(const Real &d){
        this->w/=d;
        this->s/=d;
        return *this;
    }
    Quaternion operator* (const Real &p){
        this->w *= p;
        this->s *= p;
        return *this;
    }
    Quaternion operator*(const Quaternion &p){
        float scalar = s*p.s-(w.dotProduct(p.w));
        Vec3f vector = s*p.w+p.s*w+crossProductMatrix(w)*p.w;
        return Quaternion(scalar,vector);
    }
    Quaternion operator+=(const Quaternion& a){
        s+=a.s;
        w+=a.w;
        return *this;
    }

    Mat3f rotation(){
        Real x,y,z;
        x=w.x;
        y=w.y;
        z=w.z;
        return Mat3f(
        1 - 2*y*y - 2*z*z,2*x*y - 2*s*z,2*x*z + 2*s*y,
        2*x*y + 2*s*z,1 - 2*x*x - 2*z*z,2*y*z - 2*s*x,
        2*x*z - 2*s*y,2*y*z + 2*s*x,1 - 2*x*x - 2*y*y
    );
    }
    Real normalization(){
        return sqrt(w.x*w.x+w.y*w.y+w.z*w.z+s*s);
    }
    Quaternion normalize(){
        return Quaternion(*this)/normalization();
    }
};

#endif  /* _QUATERNION_HPP_ */