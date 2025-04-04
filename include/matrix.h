#ifndef _MATRIX_
#define _MATRIX_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;
/**
 * @brief Matrix base class. The class defines a matrix of constant shape (rows and columns). Most operations will require 
 * similarily-shaped matrices and will error if the conditions are not met.
 * 
 */
class Matrix {
public:
    int n_row, n_column;
	double **data;

    // Parameterized constructor
    Matrix(const int n_row, const int n_column);
	Matrix(const int v_size);
	
	// Member operators
	double& operator () (const int row, const int column);
	double& operator () (const int n);
	Matrix& operator + (Matrix &m);
	Matrix& operator - (Matrix &m);
	Matrix& operator * (Matrix &m);
	Matrix& operator / (Matrix &m);
	Matrix& operator = (Matrix &m);

	Matrix& operator / (const double n);
	Matrix& operator * (const double n);
	Matrix& operator + (const double n);
	Matrix& operator - (const double n);

	Matrix& inv();
	Matrix& transpose();
	// Non-member operators
	friend ostream& operator << (ostream &o, Matrix &m);
};

// Operator overloading
ostream& operator << (ostream &o, Matrix &m);

// Methods
Matrix& zeros(const int n_row, const int n_column);

Matrix& zeros(const int n);

Matrix& eye(const int n_row, const int n_column);

#endif