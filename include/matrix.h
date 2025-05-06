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
	Matrix();
	
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

	Matrix& operator-();

	
	Matrix& transpose();
	Matrix& extract_vector(const int start, const int end);
	Matrix& extract_row(const int n);
	Matrix& extract_column(const int n);
	Matrix& assign_row(const int n, Matrix &m);
	Matrix& assign_column(const int n, Matrix &m);

	
	// Non-member operators
	friend ostream& operator << (ostream &o, Matrix &m);
};

// Operator overloading
ostream& operator << (ostream &o, Matrix &m);

// Methods
Matrix& zeros(const int n_row, const int n_column);

Matrix& zeros(const int n);

Matrix& union_vector(Matrix &m1, Matrix &m2);

Matrix& eye(const int n_row, const int n_column);

double dot(Matrix &m1, Matrix &m2);

double norm(Matrix &m1);

Matrix& cross(Matrix &m1, Matrix &m2);

Matrix& inv(Matrix &m);

#endif