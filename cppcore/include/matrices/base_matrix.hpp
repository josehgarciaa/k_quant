#ifndef BASEMATRIX_H
#define BASEMATRIX_H

#include <complex>
#include <vector>
#include <Eigen/Dense>
#include "typedef.hpp"



class BaseMatrix{
public:

public:



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
    virtual
    void setElement(const size_t row,const size_t  col, const std::complex<double>& value);


    // Function to set elements of the Hermitian matrix
    void addElement(int row, int col, const std::complex<double>& value) {
        if (row == col)
            data(row, col)+= value.real();
        if (row < col) 
        {
            data(row, col)+= value;
            data(col, row)+= std::conj(value);
        }
    }


    // Function to get elements of the Hermitian matrix
    std::complex<double> getElement(int row, int col) const {
        return data(row, col);
    }

    // Function to get the size (N) of the matrix
    int size() const {
        return data.rows();
    }

    // Function to perform matrix multiplication (Result is Hermitian)
    HermitianMatrix& MatScale(const HermitianMatrix& ScalMat){
        this->data = ScalMat.data * data;
        return *this;
    }

    // Function to perform matrix multiplication (Result is Hermitian)
    HermitianMatrix MatProd(const HermitianMatrix& AMat) const {
        HermitianMatrix result( this->size() );
        result.data =  AMat.data * this->data ;
        return result;
    }


    HermitianMatrix& ScaleAndSum(const double my_scale, const HermitianMatrix& A, const double a_scale){
                                this->data = my_scale*this->data + a_scale*A.data;  
                            return *this;} 

    HermitianMatrix& Square(){
                                this->data = this->data*this->data;  
                            return *this;} 


    // Function to set the matrix to the identity matrix
    HermitianMatrix& setIdentity() {
        data.setZero(); // Set all elements to zero

        // Set diagonal elements to one (identity matrix)
        for (int i = 0; i < data.rows(); ++i) {
            data(i, i) = std::complex<double>(1.0, 0.0);
        }
        return *this;
    }

    // Function to compute the trace of the matrix (Trace is real)
    double trace() const {
        return data.real().trace();
    }

    // Function to print the matrix
    void print() const {
        std::cout << data << std::endl;
    }
};



class DiagonalMatrix {
public:
    Eigen::DiagonalMatrix<complex_t, Eigen::Dynamic> diagonal;

public:
    // Constructor
    DiagonalMatrix(const size_t size = 0) : diagonal(size) {
        diagonal.diagonal().setZero();
    }

    DiagonalMatrix(std::vector< complex_t>& elems){
        if( (size_t)diagonal.diagonal().size() != (size_t)elems.size() )
            diagonal = Eigen::DiagonalMatrix<complex_t, Eigen::Dynamic>(elems.size());
        SetElements(elems);
        }


    int SetElements(const std::vector< complex_t>& elems){
        if( (size_t)diagonal.diagonal().size() != (size_t)elems.size() )
        {
            std::cerr<<"In Diagonal matrix the elements has not the same size as the allocated matrix"<<std::endl;
            exit(-1);
        }
        for(size_t idx = 0 ; idx < elems.size(); idx++)
            setElement(idx,elems[idx]);
        
        return 0;
    }


    // Function to set elements of the diagonal matrix
    inline 
    void setElement(const size_t index, const complex_t& value) {
        diagonal.diagonal()[index] = value;
    }

    // Function to add elements to the diagonal matrix
    void addElement(int index, const complex_t& value) {
        diagonal.diagonal()[index] += value;
    }

    // Function to get elements of the diagonal matrix
    std::complex<double> getElement(const size_t index) const {
        return diagonal.diagonal()[index];
    }

    // Function to get the size (N) of the matrix
    size_t size() const {
        return diagonal.diagonal().size();
    }

    // Function to perform matrix scaling (Result is diagonal)
    DiagonalMatrix& MatScale(const DiagonalMatrix& ScalMat) {
        diagonal.diagonal() = ScalMat.diagonal.diagonal().array() * diagonal.diagonal().array();
        return *this;
    }

    // Function to perform matrix multiplication (Result is diagonal)
    DiagonalMatrix MatProd(const DiagonalMatrix& AMat) const {
        DiagonalMatrix result(this->size());
        result.diagonal.diagonal() = AMat.diagonal.diagonal().array() * this->diagonal.diagonal().array();
        return result;
    }

    DiagonalMatrix& ScaleAndSum(const complex_t my_scale, const DiagonalMatrix& A, const complex_t a_scale) {
        this->diagonal.diagonal() = my_scale * this->diagonal.diagonal() + a_scale * A.diagonal.diagonal();
        return *this;
    }

    DiagonalMatrix& Square() {
        this->diagonal.diagonal() = this->diagonal.diagonal().array() * this->diagonal.diagonal().array();
        return *this;
    }

    // Function to set the matrix to the identity matrix
    DiagonalMatrix& setIdentity() {
        diagonal.diagonal().setZero();
        for (int i = 0; i < diagonal.diagonal().size(); ++i) {
            diagonal.diagonal()[i] = complex_t(1.0, 0.0);
        }
        return *this;
    }

    complex_t trace() const {
        return diagonal.diagonal().sum();
    }

    std::vector<complex_t> GetDiagonal() {
        // Create a std::vector to hold the diagonal elements
        std::vector<complex_t> diagElements(this->size());
        // Copy the diagonal elements from the Eigen matrix to the std::vector
        for (size_t i = 0; i < this->size(); ++i) {
            diagElements[i] = diagonal.diagonal()[i];
        }
    return diagElements;
    }


    // Function to print the matrix
    void print() const {
        std::cout << "Diagonal elements: " << diagonal.diagonal().transpose() << std::endl;
    }


    HermitianMatrix MatMul(const HermitianMatrix& R) const
    {
        return diagonal*R.data;
    }




};






#endif // HAMILTONIAN_READER_H
