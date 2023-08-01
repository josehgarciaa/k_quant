#include <iostream>
#include "hamiltonian.hpp"
#include "wannier90_utils.hpp"

int main() {
    // Your code using lib2 functions goes here
    std::cout << "Hello from the main program!" << std::endl;
//    auto hoppings = W90::HamiltonianFromHR("graphene_hr.dat");

    HoppingList hoppings;
    build_k_hamiltonian(0.1,0.1,0.1,hoppings);

    
    return 0;
}