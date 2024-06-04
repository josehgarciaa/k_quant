#ifndef BASE_OPERATORS_HPP
#define BASE_OPERATORS_HPP

#include <vector>
#include "typedef.hpp"

struct OperatorEntry{
    int d0,d1,d2;
    int initial_orbital; ///< Initial site of the hopping
    int final_orbital; ///< Final site of the hopping
    real_t rvalue, ivalue;
};

class Operator
{
    int basis_size;
    std::vector<OperatorEntry> list;
    
    public:
    void Reserve(const size_t size){
        list.reserve(size);
    }

    const std::vector<OperatorEntry>& GetEntryList() const 
    {
        return list;
    }    


    void AddEntry(const OperatorEntry& hop)
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

};

#endif 
