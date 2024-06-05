#ifndef BASE_OPERATORS_HPP
#define BASE_OPERATORS_HPP

#include <string>
#include <fstream>
#include <limits>
#include <iostream>
#include <limits>
#include <vector>
#include "typedef.hpp"


struct OperatorEntry{
    int d0,d1,d2;
    size_t initial_orbital; ///< Initial site of the hopping
    size_t final_orbital; ///< Final site of the hopping
    real_t rvalue, ivalue;
};

class Operator
{
    size_t basis_size;
    std::vector<OperatorEntry> list;

    protected:
    std::string EntriesToString() const
    {
        std::ostringstream oss;
        oss.precision(std::numeric_limits<double>::digits10 + 2);
        for (auto &ham_entry : this->GetEntryList()) {
            oss<<ham_entry.d0<<"    ";
            oss<<ham_entry.d1<<"    ";
            oss<<ham_entry.d2<<"    ";
            oss<<ham_entry.initial_orbital<<"    ";
            oss<<ham_entry.final_orbital<<"    ";
            oss<<ham_entry.rvalue<<"    ";
            oss<<ham_entry.ivalue<<std::endl;
        }
        return oss.str();
    }

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

    void SetBasisSize(const size_t& size)
    {
        basis_size = size;
    }

    size_t GetBasisSize() const
    {
        return basis_size;
    }

    // Method to convert the object to a string representation
    virtual std::string toString() const {
        std::ostringstream oss;
        oss.precision(std::numeric_limits<double>::digits10 + 2);
        oss << "Base Operator"<<" created using kquant"<<std::endl;
        oss << this->GetBasisSize()<<std::endl<<std::endl;
        oss << this->EntriesToString()<<std::endl;
        return oss.str();
    }

};

#endif 
