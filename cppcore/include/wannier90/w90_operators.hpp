#ifndef W90_OPERATORS_HPP
#define W90_OPERATORS_HPP

#include <vector>
#include "operators.hpp"
#include "typedef.hpp"

namespace W90{
    
    class Operator : public ::Operator {
    public:
        int num_grid_points;

        //Setters
        void SetNumGridPoints(const int& num) {  
            num_grid_points = num; }

        //Getters
        int GetNumGridPoints() const { 
            return num_grid_points; }


        


    };

    class Hamiltonian : public Operator {

    public:
        std::vector<std::string> wz_gpoints;    

        std::string GetWzGpointsLine(const int i) const 
        {
            return wz_gpoints[i];
        }

        void SetWzGpoints(const std::vector<std::string>& gpoints)
        {
            wz_gpoints = gpoints;
        }
    };
};
#endif 
