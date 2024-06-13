#ifndef W90_OPERATORS_HPP
#define W90_OPERATORS_HPP

#include <vector>
#include "operators.hpp"
#include "typedef.hpp"
#include "physics.hpp"
#include "matrices.hpp"
#include "matrices/diagonal.hpp"
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

    class Centres
    {
        size_t num_wann_centres = 0; 
        std::vector< kquant::CartVector > centres_pos;
        DiagonalMatrix u_matrix;

        public:
        Centres(){}; 

        public:
        void Reserve(const size_t size)
        {
            centres_pos.reserve(size);
        }

        size_t GetNumberOfCentres() const
        {
            return num_wann_centres;
        }
    

        void AddCenter(const kquant::CartVector& C)
        {
            centres_pos.push_back(C);
            num_wann_centres =  centres_pos.size();
        }

        void AddCentres(const std::vector< kquant::CartVector >& Cs)
        {
            centres_pos= std::vector< kquant::CartVector >(Cs);
            num_wann_centres =  centres_pos.size();
        }

        const std::vector< kquant::CartVector >& 
        GetCentres()
        {
            return centres_pos;
        }


        const kquant::CartVector& 
        GetCenter(const size_t idx) const 
        {
            return centres_pos[idx];
        }


        DiagonalMatrix Umatrix(const kquant::CartVector& k) 
        {
            kquant::CartVector x(1.,2.,3.);
            this->AddCenter(x);
            this->AddCenter(x);
            
            std::vector<complex_t> phases(num_wann_centres);
            for( size_t idx = 0 ; idx < num_wann_centres; idx++)
            {
                const auto centre = this->GetCenter(idx);
                phases[idx] = complex_t( 0, centre.dot(k));
            }
            u_matrix= DiagonalMatrix(phases);

            for(const auto& x: u_matrix.GetDiagonal() )
            {
                std::cout<<u_matrix.diagonal<<std::endl;
            }

        return u_matrix;};
    };

};

#endif 
