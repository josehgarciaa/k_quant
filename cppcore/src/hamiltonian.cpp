
#include <iostream>
#include "sparse_matrices.hpp"
#include "wannier90_utils.hpp"
#include "hamiltonian.hpp"


HermitianMatrix get_k_hamiltonian(const double k0,
                                    const double k1,
                                    const double k2,
                                    const HoppingList& hoppings)
{
    const int ham_dim = hoppings.GetBasisSize();
    HermitianMatrix k_hamiltonian(ham_dim);

    std::cout<<"ham_size"<<hoppings.GetBasisSize()<<std::endl;
    for(const auto& h: hoppings.list)
    {
            const real kr = ((real)h.d0)*k0 +
                            ((real)h.d1)*k1 +
                            ((real)h.d2)*k2 ;
            const scalar Ikr =  scalar(0,kr);
            const scalar t_expIkr =  scalar( h.rvalue,h.ivalue)*std::exp(Ikr);
            k_hamiltonian.addElement(h.initial_orbital,h.final_orbital, t_expIkr);
    }
    return k_hamiltonian;
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


