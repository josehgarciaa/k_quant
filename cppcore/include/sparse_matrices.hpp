#ifndef SPARSE_MATRICES_H
#define SPARSE_MATRICES_H

#include <complex>
#include <vector>
#include <Eigen/Dense>


typedef double real;
typedef std::complex<real> scalar;

class BlockSparse3D
{
    typedef std::vector<scalar > sparse_matrix_t;
    sparse_matrix_t data_;   
    const size_t dim0_, dim1_,dim2_, bdim_;    

    public:
    BlockSparse3D(  const size_t dim0=1,
                    const size_t dim1=1,
                    const size_t dim2=1,
                    const size_t bdim=1): 
                    dim0_(dim0),dim1_(dim1),
                    dim2_(dim2),bdim_(bdim)
    {
        const size_t size = bdim_*bdim_*dim0_*dim1_*dim2_;
        std::cout<<"The BlockSparse3D will allocate : "<<size*sizeof( std::complex<double>)/1024/1024 <<"MB of memory"<<std::endl;
        data_ = sparse_matrix_t(size);
    }

    scalar& MatrixElement(  const size_t i0,
                            const size_t i1,
                            const size_t i2,
                            const size_t i,
                            const size_t j)
    {
        return data_[  ( ( ( i2*dim1_ + i1)*dim0_ + i0)*bdim_ + i)*bdim_ + j];
    }

    scalar& MatrixBlock( const size_t i0,
                    const size_t i1,
                    const size_t i2)
    {
        return data_[  ( ( ( i2*dim1_ + i1)*dim0_ + i0)*bdim_ )*bdim_ ];
    }
};


#endif // HAMILTONIAN_READER_H
