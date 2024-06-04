#ifndef W90UTIL_PARSERS_H
#define W90UTIL_PARSERS_H

#include <string>
#include <fstream>
#include <limits>

#include "hamiltonian.hpp"

namespace W90{

/**
 * @brief Counts the total number of lines in $label_hr.dat file.
 *
 * @param filename The name of the file to be processed.
 * @return The total number of lines in the file. 
 * @throws std::runtime_error if the file cannot be opened or read.
 *
 * @note add note.
 */
unsigned long TotalLineNumber(const std::string& filename);

HoppingList HamiltonianFromHR(const std::string& filename, bool read_gridpoint);



//unsigned long HeaderLineNumber(const std::string& filename);

//unsigned long HoppingLineNumber(const std::string& filename);

};

#endif
