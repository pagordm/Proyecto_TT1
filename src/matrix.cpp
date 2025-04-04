#include "..\include\matrix.h"
/**
 * @brief Construct a new Matrix:: Matrix object
 * 
 * @param n_row The number of rows of the matrix
 * @param n_column The number of columns of the matrix
 */
Matrix::Matrix(const int n_row, const int n_column) {
    if (n_row <= 0 || n_column <= 0) {
		cout << "Matrix create: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	this->n_row = n_row;
	this->n_column = n_column;
	this->data = (double **) malloc(n_row*sizeof(double *));
	
    if (this->data == NULL) {
		cout << "Matrix create: error in data\n";
        exit(EXIT_FAILURE);
	}
	
	for(int i = 0; i < n_row; i++) {
		this->data[i] = (double *) malloc(n_column*sizeof(double));
	}
}
/**
 * @brief Accesses a element of the matrix.
 * 
 * @param row Row of the element
 * @param column Column of the element
 * @return double& A reference to the element in the row and columns
 */
double& Matrix::operator () (const int row, const int column) {
	if (row <= 0 || row > this->n_row || column <= 0 || column > this->n_column) {
		cout << "Matrix get: error in row/column\n";
        exit(EXIT_FAILURE);
	}
	
	return this->data[row - 1][column - 1];
}
/**
 * @brief Assignment operator. Copies all values of matrix m into this. Requires both matrices to have the same shape.
 * 
 * @param m The matrix to copy data from
 * @return Matrix& A reference to this, to be used in chaining assignments.
 */
Matrix& Matrix::operator = (Matrix &m) {
	if (this->n_row != m.n_row || this->n_column != m.n_column) {
		cout << "Matrix assignment: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	for(int i = 0; i < this->n_row; i++) {
		for (int j = 0; j < this->n_column; j++) {
			this->data[i][j]=m.data[i][j];
		}
	}
	return *this;
}
/**
 * @brief Add two matrices element-wise, i.e., add all values element-by-element. Both matrices must have the same shape.
 * 
 * @param m Second matrix to sum
 * @return Matrix& Result of matrix summation.
 */
Matrix& Matrix::operator + (Matrix &m) {
	if (this->n_row != m.n_row || this->n_column != m.n_column) {
		cout << "Matrix sum: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) + m(i,j);
		}
	}
	
	return *m_aux;
}
/**
 * @brief Subtracts the second matrix from the first matrix
 * 
 * @param m Matrix to subtract from this
 * @return Matrix& Resulting matrix
 */
Matrix& Matrix::operator - (Matrix &m) {
	if (this->n_row != m.n_row || this->n_column != m.n_column) {
		cout << "Matrix sub: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) - m(i,j);
		}
	}
	
	return *m_aux;
}
/**
 * @brief Prints matrix to ostream.
 * 
 * @param o ostram to print matrix to
 * @param m Matrix to print
 * @return ostream&
 */
ostream& operator << (ostream &o, Matrix &m) {
	for (int i = 1; i <= m.n_row; i++) {
        for (int j = 1; j <= m.n_column; j++)
			printf("%5.20lf ", m(i,j));
        o << "\n";
    }
	
    return o;
}
/**
 * @brief Multiplies two matrices. If matrix1.n_column!=matrix2.n_row, the method will error.
 * 
 * @param m Matrix to multiply with this
 * @return Matrix& Resulting matrix
 */
Matrix& Matrix::operator * (Matrix &m) {
	if (this->n_column != m.n_row) {
		cout << "Matrix mul: error in n_column/n_row\n";
		exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(this->n_row, m.n_column);
	
	for(int i = 1; i <= this->n_row; i++) {
		for(int j = 1; j <= m.n_column; j++) {
			(*m_aux)(i,j) = 0;
			
			for(int k = 1; k <= this->n_column; k++) {
				(*m_aux)(i,j) += (*this)(i,k) * m(k,j);
			}
		}
	}
	
	return *m_aux;
}
/**
 * @brief Multiplies this by the inverse of matrix m. Equivalent to this*m.inv()
 * 
 * @param m Matrix to invert and multiply with this.
 * @return Matrix& Resulting matrix
 */
Matrix& Matrix::operator / (Matrix &m) {
	if (this->n_column != m.n_row) {
		cout << "Matrix div: error in n_column/n_row\n";
		exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = &m.inv();
	
	return (*this) * (*m_aux);
}
/**
 * @brief Calculates the inverse of this matrix
 * 
 * @return Matrix& Reference to the inverted matrix
 */
Matrix& Matrix::inv() {
	if (this->n_row != this->n_column) {
		cout << "Matrix inv: error in n_row/n_column\n";
		exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = &eye(this->n_row, this->n_column);
	
	for(int i = 1; i <= this->n_row; i++) {
		double pivot = (*this)(i,i);

		if (pivot==0) {
			cout << "Matrix inv: singular matrix not invertible\n";
			exit(EXIT_FAILURE);
		}
		
		for(int j = 1; j <= this->n_column; j++) {
			(*this)(i,j) /= pivot;
			(*m_aux)(i,j) /= pivot;
		}
		
		for(int k = 1; k <= this->n_row; k++) {
			if (k != i) {
				double factor = (*this)(k,i);
				
				for(int j = 1; j <= this->n_column; j++) {
					(*this)(k,j) -= factor * (*this)(i,j);
					(*m_aux)(k,j) -= factor * (*m_aux)(i,j);
				}
			}
		}
	}
	
	return *m_aux;
}
/**
 * @brief Calculates the transposed matrix of this
 * 
 * @return Matrix& A reference to the transposed matrix
 */
Matrix& Matrix::transpose() {
	Matrix *m_aux = new Matrix(this->n_column, this->n_row);
	
	for(int i = 1; i <= this->n_row; i++) {
		for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(j,i) = (*this)(i,j);
		}
	}
	
	return *m_aux;
}
/**
 * @brief Generates a matrix with all zeros.
 * 
 * @param n_row Rows of the matrix
 * @param n_column Columns of the matrix
 * @return Matrix& A reference to the matrix
 */
Matrix& zeros(const int n_row, const int n_column) {
	Matrix *m_aux = new Matrix(n_row, n_column);
	
	for(int i = 1; i <= n_row; i++) {
		for(int j = 1; j <= n_column; j++) {
			(*m_aux)(i,j) = 0;
		}
	}
	
	return (*m_aux);
}
/**
 * @brief Generates a matrix with all ones.
 * 
 * @param n_row Rows of the matrix
 * @param n_column Columns of the matrix
 * @return Matrix& A reference to the matrix
 */
Matrix& eye(const int n_row, const int n_column) {
	Matrix *m_aux = new Matrix(n_row, n_column);
	
	for(int i = 1; i <= n_row; i++) {
		for(int j = 1; j <= n_column; j++) {
			if (i == j)
				(*m_aux)(i,j) = 1;
			else
				(*m_aux)(i,j) = 0;
		}
	}
	
	return (*m_aux);
}
/**
 * @brief Divides a matrix element-wise by a double
 * 
 * @param n Double to divide the matrix by.
 * @return Matrix& Reference to the resulting matrix.
 */
Matrix& Matrix::operator / (const double n) {
	
	if (n==0) {
		cout << "Double inv: error in n\n";
		exit(EXIT_FAILURE);
	}
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);

	for(int i = 1; i <= this->n_row; i++) {
		for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j)/n;
		}
	}
	return *m_aux;

}

/**
 * @brief Multiplies a matrix element-wise by a double
 * 
 * @param n Double to multiply the matrix by.
 * @return Matrix& Reference to the resulting matrix.
 */
Matrix& Matrix::operator * (const double n) {
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);

	for(int i = 1; i <= this->n_row; i++) {
		for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j)*n;
		}
	}
	return *m_aux;

}
/**
 * @brief Adds a matrix element-wise by a double
 * 
 * @param n Double to add the matrix by.
 * @return Matrix& Reference to the resulting matrix.
 */
Matrix& Matrix::operator + (const double n) {
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);

	for(int i = 1; i <= this->n_row; i++) {
		for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j)+n;
		}
	}
	return *m_aux;

}