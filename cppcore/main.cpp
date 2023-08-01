#include <iostream>
#include <vector>
#include "sparse_matrices.hpp"
#include "hamiltonian.hpp"
#include "wannier90_utils.hpp"

int main() {
    // Your code using lib2 functions goes here
    std::cout << "Hello from the main program!" << std::endl;
    auto hoppings = W90::HamiltonianFromHR("graphene_hr.dat");
    
    const auto Ham = get_k_hamiltonian(0.1,0.1,0.1,hoppings);
    const int ham_size = Ham.size();

    #pragma omp parallel for
    for(int nkpt = 0; nkpt < 10000; nkpt++ )
    {
        // Chebyshev module
        const int M = 100;
        std::vector< HermitianMatrix > ChebT(M);
        for(int m=0; m < M; m ++)
            ChebT[m] =  HermitianMatrix(ham_size).setIdentity();

        ChebT[1] = Ham;
        for( int m = 2; m < M; m++)
        { 
            ChebT[m] = ChebT[m-1].MatProd(Ham);
            ChebT[m].ScaleAndSum(2, ChebT[m-2],-1);
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


    
    return 0;
}