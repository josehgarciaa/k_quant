#ifndef DIAGONALMATRIX_H
#define DIAGONALMATRIX_H

#include <complex>
#include <vector>
#include "matrices/base_matrix.hpp"


class DiagonalMatrix : public BaseMatrix{
public:
    Eigen::DiagonalMatrix<complex_t, Eigen::Dynamic> diagonal;

public:
    // Constructor
    DiagonalMatrix(const size_t size = 0) : DiagonalMatrix(),diagonal(size) {
        this->setMatrixType("Diagonal Matrix");
        diagonal.diagonal().setZero();}

    DiagonalMatrix(std::vector< complex_t>& elems){
        this->setMatrixType("Diagonal Matrix");
        if( (size_t)diagonal.diagonal().size() != (size_t)elems.size() )
            diagonal = Eigen::DiagonalMatrix<complex_t, Eigen::Dynamic>(elems.size());
        this->SetElements(elems);}


/*************************BASE CLASS FUNCTIONS*********************/
    

    //<<<<<<<<<<<<<<<<<<<GETTERS>>>>>>>>>>>>>>>>>>>>>>>>>//
    const std::complex<double>& 
    getElement(size_t row, size_t col) const {
        if( row == col )
            return diagonal.diagonal()[row];     
        else
            return this->null_complex;}

    size_t size() const {
        return diagonal.diagonal().size();}

    //<<<<<<<<<<<<<<<<<<<SETTERS>>>>>>>>>>>>>>>>>>>>>>>>>//
    BaseMatrix&  
    SetElements(const std::vector<std::complex<double>>& values){
        if( (size_t)diagonal.diagonal().size() != (size_t)elems.size() )
        {
            std::cerr<<"In Diagonal matrix the elements has not the same size as the allocated matrix"<<std::endl;
            exit(-1);
        }
        for(size_t idx = 0 ; idx < elems.size(); idx++)
            setElement(idx,elems[idx]);       
        return *this;}

    BaseMatrix& 
    setElement(const size_t row,const size_t col, const complex_t& value){ 
        if (row == col)
            diagonal.diagonal()[row] = value; 
        return this*}

    //<<<<<<<<<<<<<<<<<<<MATRIX OPERATIONS>>>>>>>>>>>>>>>>>>>>>>>>>//

    // Function to compute the trace of the matrix (Trace is real)
    std::complex<double> trace() const{
        return diagonal.diagonal().size();}

    //<<<<<<<<<<<<<<<<<<<IO OPERATIONS>>>>>>>>>>>>>>>>>>>>>>>>>//

    // Function to print the matrix
    void print() const{
        std::cerr<<"The print method is not define in "<<this->GetMatrixType()<<" class"<<std::endl;}

    /*************************DERIVE CLASS FUNCTIONS*********************/

    // Function to set elements of the diagonal matrix
    inline 
    DiagonalMatrix& setElement(const size_t index, const complex_t& value){ 
        diagonal.diagonal()[index] = value; 
        return *this;}

    // Function to get elements of the diagonal matrix
    const std::complex<double>& 
    getElement(const size_t index) const{ 
        return diagonal.diagonal()[index];}

    // Function to perform matrix scaling (Result is diagonal)
    DiagonalMatrix& MatScale(const DiagonalMatrix& ScalMat){
        diagonal.diagonal() = ScalMat.diagonal.diagonal().array() * diagonal.diagonal().array();
        return *this;}

    // Function to perform matrix multiplication (Result is diagonal)
    DiagonalMatrix& MatProd(const DiagonalMatrix& AMat) const{
        DiagonalMatrix result(this->size());
        result.diagonal.diagonal() = AMat.diagonal.diagonal().array() * this->diagonal.diagonal().array();
        return result;}

    DiagonalMatrix& ScaleAndSum(const complex_t my_scale, const DiagonalMatrix& A, const complex_t a_scale) {
        this->diagonal.diagonal() = my_scale * this->diagonal.diagonal() + a_scale * A.diagonal.diagonal();
        return *this;}

    DiagonalMatrix& Square(){
        this->diagonal.diagonal() = this->diagonal.diagonal().array() * this->diagonal.diagonal().array();
        return *this;}

    // Function to set the matrix to the identity matrix
    DiagonalMatrix& setIdentity(){
        diagonal.diagonal().setZero();
        for (int i = 0; i < diagonal.diagonal().size(); ++i) {
            diagonal.diagonal()[i] = complex_t(1.0, 0.0);
        }
        return *this;}

    std::vector<complex_t> GetDiagonal(){
        // Create a std::vector to hold the diagonal elements
        std::vector<complex_t> diagElements(this->size());
        // Copy the diagonal elements from the Eigen matrix to the std::vector
        for (size_t i = 0; i < this->size(); ++i) {
            diagElements[i] = diagonal.diagonal()[i];
        }
    return diagElements;
    }

};

#endif 