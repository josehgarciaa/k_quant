#ifndef BASEMATRIX_H
#define BASEMATRIX_H

#include <complex>
#include <vector>
#include <Eigen/Dense>
#include "typedef.hpp"

class BaseMatrix{

    /**
     * @brief Set the value of an element in the matrixlocated at ith row and jth column using base-zero indexing.
     *
     * Set the value of an element in the matrix located at 
     * the ith row and jth column using base-zero indexing.
     * 
     *
     * @param[in] i Description of the first argument (e.g., an integer).
     * @param[in] j Description of the second argument (e.g., a double).
     * @param[in] value Description of the third argument (e.g., a string).
     * @return Pointer to the object 
     *
     * @note Additional notes about the function, if any.
     * @warning Any warnings about the function usage.
     * @pre Precondition that must be satisfied before calling this function.
     * @post Postcondition guaranteed after the function returns.
     * @exception Description of any exceptions thrown by the function.
     * @see Reference to other related functions or documentation.
     */

    //<<<<<<<<<<<<<<<<<<<CONSTRUCTORS>>>>>>>>>>>>>>>>>>>>>>>>>//

    public:
    BaseMatrix(const std::string& name){
        matrix_type_="BaseMatrix";}


    //<<<<<<<<<<<<<<<<<<<GETTERS>>>>>>>>>>>>>>>>>>>>>>>>>//

    // Function to get elements of the Hermitian matrix
    virtual
    const std::complex<double>& 
    getElement(const size_t row,const size_t col) const {
        std::cerr<<"The getElement method is not define in "<<this->GetMatrixType()<<" class"
                 <<"returning: "<<this->null_complex
                 <<std::endl;
        return this->null_complex;}

    // Function to get the size (N) of the matrix
    virtual
    size_t size() const {
        std::cerr<<"The size method is not define in "
                 <<this->GetMatrixType()<<" class"
                 <<"returning: "<<this->null_size
                 <<std::endl;
        return this->null_size;}

    //<<<<<<<<<<<<<<<<<<<SETTERS>>>>>>>>>>>>>>>>>>>>>>>>>//
    virtual
    BaseMatrix& setElements(const std::vector<std::complex<double>>& values){
        std::cerr<<"The setElements method is not define in "
                 <<this->GetMatrixType()<<" class"
                 <<std::endl;  
        return *this;}

    virtual
    BaseMatrix& setElement(const size_t row,const size_t  col, const std::complex<double>& value){
        std::cerr<<"The setElement method is not define in "
                 <<this->GetMatrixType()<<" class"
                 <<std::endl;
        return *this;}

    //<<<<<<<<<<<<<<<<<<<MATRIX OPERATIONS>>>>>>>>>>>>>>>>>>>>>>>>>//

    // Function to compute the trace of the matrix (Trace is real)
    virtual
    std::complex<double> trace() const{
        std::cerr<<"The trace method is not define in "<<this->GetMatrixType()
                 <<" class. Returning"<<this->null_complex
                 <<std::endl;
        return this->null_complex;}


    //<<<<<<<<<<<<<<<<<<<IO OPERATIONS>>>>>>>>>>>>>>>>>>>>>>>>>//

    // Function to print the matrix
    virtual
    void print() const{
        std::cerr<<"The print method is not define in "<<this->GetMatrixType()<<" class"<<std::endl;}

    //Protected & Privates
    protected:
        std::string matrix_type_;
        double null_double=0.0;
        size_t null_size=0;
        std::complex<double> null_complex = {0.,0.};
        std::vector<std::complex<double>> null_vecCD;

        inline void setMatrixType(const std::string matrix_type) {matrix_type_ = matrix_type;}
        inline const std::string& GetMatrixType() const {return matrix_type_;}

};







#endif 