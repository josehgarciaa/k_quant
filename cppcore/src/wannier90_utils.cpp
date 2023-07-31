#include "wannier90_utils.hpp"



unsigned long W90::TotalLineNumber(const std::string& filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error( "Cannot open the file: "+filename );
  }
  unsigned long line_count = 0;
  while (file.ignore(std::numeric_limits<std::streamsize>::max(), '\n')) {
    ++line_count;
  }
  file.close();
  return line_count;
}


HoppingList W90::HamiltonianFromHR(const std::string& filename) {
    HoppingList hops;
    
    //Reserve enough continious space for the hoppings
    //by using the total number of lines of the file
    long int line_number = W90::TotalLineNumber(filename);
    hops.list.reserve(line_number);
    
    //Required to map the string into numbers
    std::istringstream iss;

    // Open the file and read each string with a fixed precision
    std::ifstream input_file(filename);
    if (!input_file.is_open()) {
        throw std::runtime_error("Cannot open the file: "+filename);
    }
    input_file.precision(std::numeric_limits<double>::digits10 + 2);

    // Read the wannier file following wannier90 formatting
    Hopping hop; //variable to storage the hopping     
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
                hops.AddHopping(hop);
            }    
        }
    }
    input_file.close();

    return hops;
}


