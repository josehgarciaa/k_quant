#ifndef MATRICES_H
#define MATRICES_H

#include <complex>
#include <vector>
#include <Eigen/Dense>



class HermitianMatrix {
private:
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

class DiagonalMatrix {
private:
    Eigen::DiagonalMatrix<std::complex<double>, Eigen::Dynamic> diagonal;

public:
    // Constructor
    DiagonalMatrix(const int size = 0) : diagonal(size) {
        diagonal.diagonal().setZero();
    }

    // Function to set elements of the diagonal matrix
    void setElement(int index, const std::complex<double>& value) {
        diagonal.diagonal()[index] = value;
    }

    // Function to add elements to the diagonal matrix
    void addElement(int index, const std::complex<double>& value) {
        diagonal.diagonal()[index] += value;
    }

    // Function to get elements of the diagonal matrix
    std::complex<double> getElement(int index) const {
        return diagonal.diagonal()[index];
    }

    // Function to get the size (N) of the matrix
    int size() const {
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

    DiagonalMatrix& ScaleAndSum(const double my_scale, const DiagonalMatrix& A, const double a_scale) {
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
            diagonal.diagonal()[i] = std::complex<double>(1.0, 0.0);
        }
        return *this;
    }

    // Function to compute the trace of the matrix (Trace is real)
    double trace() const {
        return diagonal.diagonal().real().sum();
    }

    // Function to print the matrix
    void print() const {
        std::cout << "Diagonal elements: " << diagonal.diagonal().transpose() << std::endl;
    }
};






#endif // HAMILTONIAN_READER_H
