#ifndef HAMILTONIAN_READER_H
#define HAMILTONIAN_READER_H

#include <fstream>
#include <iostream>
#include <complex>
#include <vector>
#include <array>
#include<limits>
#include "matrices.hpp"
#include "sparse_matrices.hpp"
#include "typedef.hpp"


struct OperatorEntry{
    int d0,d1,d2;
    int initial_orbital; ///< Initial site of the hopping
    int final_orbital; ///< Final site of the hopping
    real_t rvalue, ivalue;
};

struct Operator
{
    int basis_size;
    int num_grid_points;
    std::vector<std::string> wz_gpoints;    
    std::vector<Hopping> list;
    
    void AddHopping(const Hopping& hop)
    {
        list.push_back(hop);
    }

    void SetBasisSize(const int& size)
    {
        basis_size = size;
    }

    int GetBasisSize() const
    {
        return basis_size;
    }

    void SetNumGridPoints(const int& num)
    {
        num_grid_points = num;
    }

    int GetNumGridPoints() const
    {
        return num_grid_points;
    }

    void SetWzGpoints(const std::vector<std::string>& gpoints)
    {
        wz_gpoints = gpoints;
    }

    std::vector<std::string>& GetWzGpoints() const 
    {
        return wz_gpoints;
    }


    // Conversion methods
    std::string toString() const {
        std::ostringstream oss;
        oss.precision(std::numeric_limits<double>::digits10 + 2);
        oss << "Hamiltonian created using kquant"<<std::endl;
        oss << this->GetBasisSize()<<std::endl;
        oss << this->GetNumGridPoints()<<std::endl;
        for (auto &line : this->GetWzGpoints()) 
            oss << line<<std::endl;
//        oss << this->EntriesToString()<<std::endl;

        return oss.str();
    }


};

/**
 * @brief Structure for an entry in the Hamiltonian.
 */
struct HamiltonianEntry {
    int initial_site; ///< Initial site of the hopping
    int final_site; ///< Final site of the hopping
    std::complex<double> hopping_energy; ///< Hopping energy, as a complex number
    std::vector<int> displacement_vector;  ///< Displacement vector, assuming 3D vector
};

/**
 * @brief Reads the Hamiltonian from a file.
 * @param filename The name of the file to read from.
 * @return A vector of HamiltonianEntry, representing the Hamiltonian.
 *
 * This function opens the file named filename, and reads the Hamiltonian data into a vector of HamiltonianEntry structures.
 * Each line of the file should represent a single entry in the Hamiltonian and should be formatted as follows:
 * <initial_site> <final_site> <hopping_real_part> <hopping_imaginary_part> <displacement_vector>
 */
std::vector<HamiltonianEntry> readHamiltonian(const std::string& filename);

unsigned long CountLines(const std::string& filename);


HermitianMatrix get_onsites(    const HoppingList& hoppings,
                                const HoppingList& positions);


HermitianMatrix get_k_hamiltonian(  const double k0,
                                    const double k1,
                                    const double k2,
                                    const HoppingList& hoppings,
                                    const HoppingList& positions);


#endif // HAMILTONIAN_READER_H
