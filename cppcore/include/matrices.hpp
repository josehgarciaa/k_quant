#ifndef MATRICES_H
#define MATRICES_H

#include <complex>
#include <vector>
#include <Eigen/Dense>
#include "typedef.hpp"

#include "matrices/base_matrix.hpp"


class HermitianMatrix {
public:
    Eigen::MatrixXcd data;

public:
    // Constructor
    HermitianMatrix(const int size=0) : data(size, size) {        // Initialize the matrix to zero
        data.setZero();}

    // Function to set elements of the Hermitian matrix
    void setElement(int row, int col, const std::complex<double>& value) {
        if (row == col)
            data(row, col) = value.real();
        if (row < col) 
        {
            data(row, col) = value;
            data(col, row) = std::conj(value);
        }
    }

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










#endif // HAMILTONIAN_READER_H
