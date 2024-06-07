#ifndef W90_PARSERS_H
#define W90_PARSERS_H

#include <string>
#include <fstream>
#include <limits>
#include <iostream>

#include "typedef.hpp"
#include "physics.hpp"
#include "operators.hpp"
#include "w90_operators.hpp"



namespace W90{

    /**
     * @brief Counts the total number of lines in $label_hr.dat file.
     *
     * @param filename The name of the file to be processed.
     * @return The total number of lines in the file. 
     * @throws std::runtime_error if the file cannot be opened or read.
     *
     * @note add note.
     **/
    size_t TotalLineNumber(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error( "Cannot open the file: "+filename );
        }
        size_t line_count = 0;
        while (file.ignore(std::numeric_limits<std::streamsize>::max(), '\n')) {
            ++line_count;
        }
        file.close();
        return line_count;
    };

    /**
     * @brief Counts the total number of lines in $label_hr.dat file.
     *
     * @param filename The name of the file to be processed.
     * @return The total number of lines in the file. 
     * @throws std::runtime_error if the file cannot be opened or read.
     *
     * @note add note.
     */
    W90::Hamiltonian ReadHamiltonian(const std::string& filename)
    {
        W90::Hamiltonian hops;
        
        //Reserve enough continious space for the hoppings
        //by using the total number of lines of the file
        size_t line_number = W90::TotalLineNumber(filename);
        hops.Reserve(line_number);
    
        //Required to map the string into numbers
        std::istringstream iss;

        // Open the file and read each string with a fixed precision
        std::ifstream input_file(filename);
        if (!input_file.is_open()) {
            throw std::runtime_error("Cannot open the file: "+filename);
        }
        input_file.precision(std::numeric_limits<double>::digits10 + 2);

        // Read the wannier file following wannier90 formatting
        OperatorEntry hop; //variable to storage the hopping     
        bool read_basis_size=false;
        bool read_grid_points=false;
        bool read_wz_gpoints=false;

        std::string line,comment;
        std::getline(input_file, comment); //The first line is a comment
        while (std::getline(input_file, line)) {
            if (!line.empty() && line[0] != '#'){
                if( !read_basis_size )
                {
                    try {
                        hops.SetBasisSize( std::stoi(line) );
                        } 
                    catch(std::exception const & e)
                    {
                        std::cerr   << "Error while converting line:"<<line<<"from file" << filename<<" into basis size"<<std::endl;
                        exit(-1);
                    }
                    read_basis_size = true;
                }
                else if( !read_grid_points )
                {
                    try {
                        hops.SetNumGridPoints( std::stoi(line) );
                        } 
                    catch(std::exception const & e)
                    {
                        std::cerr << "Error while reading the number of grid points" << e.what() 
                                << "from file" << filename<< std::endl;
                        exit(-1);
                    }
                    read_grid_points = true;
                }
                else  if (!read_wz_gpoints)
                {
                    std::vector<std::string> wz_gpoints;
                    // By convention there is always 15 grid points per line
                    const int num_grid_points_lines = std::ceil( (double)hops.GetNumGridPoints() / 15.0);
                    for (int n = 0; n < num_grid_points_lines; ++n) {
                    wz_gpoints.push_back(line);
                    std::getline(input_file, line);
                    }
                    hops.SetWzGpoints(wz_gpoints);
                    read_wz_gpoints = true;
                }
                else 
                {
                    iss.str(line);
                    iss.clear();   
                    iss >> hop.d0 
                        >> hop.d1 
                        >> hop.d2 
                        >> hop.initial_orbital 
                        >> hop.final_orbital
                        >> hop.rvalue >> hop.ivalue;
                    hop.initial_orbital-=1; 
                    hop.final_orbital-=1; 
                    hops.AddEntry(hop);
                }    
            }
        }
        input_file.close();
        return hops;
    };

    /**
     * @brief Counts the total number of lines in $label_hr.dat file.
     *
     * @param filename The name of the file to be processed.
     * @return The total number of lines in the file. 
     * @throws std::runtime_error if the file cannot be opened or read.
     *
     * @note add note.
     */
    std::vector<W90::Position> RedPositionOperator(const std::string& filename){
        std::vector<W90::Position> entries_pos(3);
        
        //Reserve enough continious space for the hoppings
        //by using the total number of lines of the file
        size_t line_number = W90::TotalLineNumber(filename);
        for(size_t i=0; i < CartDIM; i++)
            entries_pos[i].Reserve(line_number);
        
        //Required to map the string into numbers
        std::istringstream iss;

        // Open the file and read each string with a fixed precision
        std::ifstream input_file(filename);
        if (!input_file.is_open()) {
            throw std::runtime_error("Cannot open the file: "+filename);
        }
        input_file.precision(std::numeric_limits<double>::digits10 + 2);

        // Read the wannier file following wannier90 formatting
        OperatorEntry entry; //variable to storage the hopping     
        bool read_basis_size=false;
        bool read_grid_points=false;

        std::string line,comment;
        std::getline(input_file, comment); //The first line is a comment
        while (std::getline(input_file, line)) {
            if (!line.empty() && line[0] != '#'){
                if( !read_basis_size )
                {
                    try {    
                            for(size_t i=0; i < CartDIM; i++)
                                entries_pos[i].SetBasisSize( std::stoi(line) );
                        } 
                    catch(std::exception const & e)
                    {
                        std::cerr   << "Error while converting line:"<<line<<"from file" << filename<<" into basis size"<<std::endl;
                        exit(-1);
                    }
                    read_basis_size = true;
                }
                else if( !read_grid_points )
                {
                    try {
                            for(size_t i=0; i < CartDIM; i++)
                                entries_pos[i].SetNumGridPoints( std::stoi(line) );
                        } 
                    catch(std::exception const & e)
                    {
                        std::cerr << "Error while reading the number of grid points" << e.what() 
                                << "from file" << filename<< std::endl;
                        exit(-1);
                    }
                    read_grid_points = true;
                }
                else 
                {
                    iss.str(line);
                    iss.clear();
                    real_t rvalue, ivalue;   
                    iss >> entry.d0 
                        >> entry.d1 
                        >> entry.d2 
                        >> entry.initial_orbital 
                        >> entry.final_orbital;
                    entry.initial_orbital-=1; 
                    entry.final_orbital-=1; 
                    for(size_t i=0; i < CartDIM; i++)
                    {
                        //We use this loop to read the imaginary part for each X,Y,Z cartesian direction
                        iss >> rvalue >> ivalue;
                        entry.rvalue =rvalue;
                        entry.ivalue =ivalue;
                        entries_pos[i].AddEntry(entry);
                    }                
                }    
            }
        }
        input_file.close();
        return entries_pos;
    };


    int WriteW90CentersFrom(const std::vector<W90::Position>& PosOp, const std::string filename){
        std::ofstream ofs(filename);
        std::vector< std::vector<real_t> > centres(CartDIM);
        for(size_t i=0; i < CartDIM; i++)
        {
            centres[i] = PosOp[i].GetOnsitePositions();
        }
        if( PosOp[0].GetBasisSize()!= PosOp[1].GetBasisSize() ||
            PosOp[1].GetBasisSize()!= PosOp[2].GetBasisSize()){
            std::cerr<<"In PositionOpToW90Centers function the size of the Position operator is not the same "<<std::endl;
            }

        ofs<<"    "<<PosOp[0].GetBasisSize()<<std::endl;
        ofs<<" Wannier centres, written by kquant"<<std::endl;
        for(size_t n =0; n < PosOp[0].GetBasisSize(); n++ )
        {
            ofs<<"X"<<"\t "<<centres[0][n]<<"\t "<<centres[1][n]<<"\t "<<centres[2][n]<<"\t "<<std::endl;
        }
        ofs.close();
    return 0;}


    std::vector< std::vector<real_t> >  ReadW90CentersFrom(const std::string filename)
    {
        std::vector< std::vector<real_t> > centres;
        //Required to map the string into numbers
        std::istringstream iss;
        // Open the file and read each string with a fixed precision
        std::ifstream input_file(filename);
        if (!input_file.is_open()) {
            throw std::runtime_error("Cannot open the file: "+filename);
        }
        input_file.precision(std::numeric_limits<double>::digits10 + 2);

        std::string line,comment;
        bool read_num_sites = false;
        bool read_comment = false;        
        while (std::getline(input_file, line)) {
            if (!line.empty() && line[0] != '#'){
                if( !read_num_sites )
                {
                    std::cout<<line<<std::endl;
                    try {    
                            const size_t num_sites = (size_t)std::stoul(line);
                            centres.reserve(num_sites);                            
                        } 
                    catch(std::exception const & e)
                    {
                        std::cerr   << "Error while converting line:"<<line<<"from file" << filename<<" into num_sites"<<std::endl;
                        exit(-1);
                    }
                    read_num_sites = true;
                }
                else if( !read_comment )
                {
                    comment = line;
                    read_comment = true;
                }
                else 
                {
                    iss.str(line);
                    iss.clear();
                    std::string X;
                    std::vector<real_t> R(3);   
                    iss >> X 
                        >> R[0]
                        >> R[1] 
                        >> R[2];
                    if ( X == "X" ) 
                    {
                        centres.push_back(R);
                    }
                }    
            }
        }
        input_file.close();
        return centres;}
        
};



#endif // HAMILTONIAN_READER_H
