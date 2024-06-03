#include <iostream>
#include <vector>
#include "wannier90_utils.hpp"
#include "hamiltonian.hpp"


int main(int argc, char* argv[]) {
    if (argc != 2 )
    {
        std::cout << "Please pass the label of the file to read" << std::endl;
        return -1;
    }
    else
        std::cout << "Reading: " << argv[1] <<std::endl;

    auto syst_name = std::string(argv[1]);
    HoppingList ham_entries = W90::HamiltonianFromHR(syst_name+"_hr.dat",true); 
    HoppingList pos_entries = W90::HamiltonianFromHR(syst_name+"_r.dat",false); 
    

    int scdim[3];
    scdim[0]= 1;
    scdim[1]= 1;
    scdim[2]= 1;

    for( int k0=0;k0<scdim[0]; k0++)
        for( int k1=0;k1<scdim[1]; k1++)
            for( int k2=0;k2<scdim[2]; k2++)
            {
                const auto Ham = get_k_hamiltonian(k0,k1,k2,ham_entries, pos_entries);
            }
    
//    std::cout<<ham_size<<std::endl<<std::endl;


    return 0;
}