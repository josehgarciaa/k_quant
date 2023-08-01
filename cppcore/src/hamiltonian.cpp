
#include <iostream>
#include "sparse_matrices.hpp"
#include "wannier90_utils.hpp"
#include "hamiltonian.hpp"



void build_k_hamiltonian(const double kx,
                         const double ky,
                         const double kz,
                         const HoppingList& hoppings)
{

/*
    // Fill the matrices (Replace these values with your actual matrix data)
    A << std::complex<double>(1.0, 2.0), std::complex<double>(3.0, 4.0), std::complex<double>(5.0, 6.0),
         std::complex<double>(7.0, 8.0), std::complex<double>(9.0, 10.0), std::complex<double>(11.0, 12.0),
         std::complex<double>(13.0, 14.0), std::complex<double>(15.0, 16.0), std::complex<double>(17.0, 18.0);

    B << std::complex<double>(2.0, 1.0), std::complex<double>(4.0, 3.0), std::complex<double>(6.0, 5.0),
         std::complex<double>(8.0, 7.0), std::complex<double>(10.0, 9.0), std::complex<double>(12.0, 11.0),
         std::complex<double>(14.0, 13.0), std::complex<double>(16.0, 15.0), std::complex<double>(18.0, 17.0);
*/

/*
    #pragma omp parallel for
    for(size_t n=0; n< kpoint_number; n++)
    {
        const size_t i0= n%kdim0;
        const size_t i1= (n/kdim0)%kdim1;
        const size_t i2= (n/kdim0/kdim1)%kdim2;
        for(const auto& h: hoppings.list)
        {
            const real k0 = (real)i0/(real)kdim0;
            const real k1 = (real)i1/(real)kdim1;
            const real k2 = (real)i2/(real)kdim2;
            
            const real kr = ((real)h.d0)*k0 +
                            ((real)h.d1)*k1 +
                            ((real)h.d2)*k2 ;
            const scalar Ikr =  scalar(0,kr);
            const scalar t_expIkr =  scalar( h.rvalue,h.ivalue)*std::exp(Ikr);
            hamiltonian_matrix.MatrixElement(i0,i1,i2,h.site_i,h.site_j) += t_expIkr;
        }
    }
*/
}





void build_hamiltonian(std::string system_label)
{
    auto hoppings = W90::HamiltonianFromHR(system_label);
    const size_t kdim0 = 100;
    const size_t kdim1 = 100;
    const size_t kdim2 = 1;
    const size_t odim = hoppings.basis_size;
    const size_t size = odim*odim*kdim0*kdim1*kdim2;

    auto hamiltonian_matrix = BlockSparse3D(kdim0,kdim1,kdim2, odim);

    const size_t kpoint_number = kdim0*kdim1*kdim2;
    #pragma omp parallel for
    for(size_t n=0; n< kpoint_number; n++)
    {
        const size_t i0= n%kdim0;
        const size_t i1= (n/kdim0)%kdim1;
        const size_t i2= (n/kdim0/kdim1)%kdim2;
        for(const auto& h: hoppings.list)
        {
            const real k0 = (real)i0/(real)kdim0;
            const real k1 = (real)i1/(real)kdim1;
            const real k2 = (real)i2/(real)kdim2;
            
            const real kr = ((real)h.d0)*k0 +
                            ((real)h.d1)*k1 +
                            ((real)h.d2)*k2 ;
            const scalar Ikr =  scalar(0,kr);
            const scalar t_expIkr =  scalar( h.rvalue,h.ivalue)*std::exp(Ikr);
            hamiltonian_matrix.MatrixElement(i0,i1,i2,h.initial_orbital,h.final_orbital) += t_expIkr;
        }
    }
}


    // Create complex matrices A and B
//     HermitianMatrix Tm0(N), Tm1(N), Tm2(N);
//    for(int i=0; i <N; i++) 
//    for(int j=0; j <N; j++)
//    { 
//        Tm0.setElement(i, j, std::complex<double>(i,j) );
//        Tm1.setElement(i, j, std::complex<double>(1,0) );
//        Tm2.setElement(i, j, std::complex<double>(0,1) );
//    }
 //   Tm0.ScaleAndSum(1, Tm1, 1).print();
//    Tm0.Square().print();
