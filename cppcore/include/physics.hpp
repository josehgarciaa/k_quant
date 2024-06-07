#ifndef PHYSICS_H
#define PHYSICS_H

#include <Eigen/Dense>

enum class Cart{ 
    X, Y, Z, DIM
};

const size_t CartDIM =3;


namespace kquant{

using namespace Eigen;

class CartVector: public Vector3d {

    public:
    CartVector(double x, double y, double z) : Vector3d(x,y,z)
    {
        data3d = Eigen::Vector3d(x,y,z);
    }

    inline
    double dot(const CartVector& xr) const
    {
        return data3d.adjoint() *(xr.data3d);   
    }



    Eigen::Vector3d data3d;



};



};

#endif // HAMILTONIAN_READER_H



