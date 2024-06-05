#ifndef W90_OPERATORS_HPP
#define W90_OPERATORS_HPP

#include <vector>
#include "operators.hpp"
#include "typedef.hpp"
#include <iostream>

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

        // Conversion methods
        std::string toString() const {
            std::ostringstream oss;
            oss.precision(std::numeric_limits<double>::digits10 + 2);
            oss << "Wannier90 Operator created using kquant"<<std::endl;
            oss << this->GetBasisSize()<<std::endl;
            oss << this->GetNumGridPoints()<<std::endl<<std::endl;
            oss << this->EntriesToString()<<std::endl;
            return oss.str();
        };
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


        const std::vector<std::string>& GetWzGpoints() const 
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
            oss << this->EntriesToString()<<std::endl;
            return oss.str();
        }

    };

    class Position : public Operator {

    public:
        
        std::vector<real_t> GetOnsitePositions() const
        {
            const double num_pos = this->GetBasisSize();
            std::vector<real_t> positions(num_pos);

            int onsite_count = 0;
            for (auto &entry : this->GetEntryList())
                if( entry.d0==0 && entry.d1==0&& entry.d2==0&&
                    entry.initial_orbital == entry.final_orbital){
                    const size_t p = entry.initial_orbital; 
                    positions[p]= entry.rvalue;
                    onsite_count++;
                    };
            if (onsite_count != num_pos)
            {
                std::cerr<<"The number of onsites and positions does not match"<<std::endl;
            }
        return positions;};


        // Conversion methods
        std::string toString() const {
            std::ostringstream oss;
            oss.precision(std::numeric_limits<double>::digits10 + 2);
            oss << "Position created using kquant"<<std::endl;
            oss << this->GetBasisSize()<<std::endl;
            oss << this->EntriesToString()<<std::endl;
        return oss.str();};

    };

};

#endif 
