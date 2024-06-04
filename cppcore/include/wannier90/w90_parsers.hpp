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
    int WriteHamiltonian(const W90::Hamiltonian& hops, const std::string filename)
    {
        /*std::cout<<"Created using k-quant"<<std::endl;
        std::cout<<"          "<<hops.GetBasisSize()<<std::endl;
        std::cout<<"          "<<hops.GetNumGridPoints()<<std::endl;
        const int num_grid_points_lines = std::ceil( (double)hops.GetNumGridPoints() / 15.0);
        for (int i=0; i <num_grid_points_lines; i++ )
            std::cout<< hops.GetWzGpointsLine(i)<<std::endl;

        for (auto &ham_entry : hops.GetEntryList()) {
            std::cout<<ham_entry.d0<<"    ";
            std::cout<<ham_entry.d1<<"    ";
            std::cout<<ham_entry.d2<<"    ";
            std::cout<<ham_entry.initial_orbital<<"    ";
            std::cout<<ham_entry.final_orbital<<"    ";
            std::cout<<ham_entry.rvalue<<"    ";
            std::cout<<ham_entry.ivalue<<std::endl;
        }*/
        return 0;
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
    std::vector<W90::Operator> RedPositionOperator(const std::string& filename){
        std::vector<W90::Operator> entries_pos(3);
        
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
};



#endif // HAMILTONIAN_READER_H
