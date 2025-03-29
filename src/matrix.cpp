#include "..\include\matrix.h"

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

ostream& operator << (ostream &o, Matrix &m) {
	for (int i = 1; i <= m.n_row; i++) {
        for (int j = 1; j <= m.n_column; j++)
			printf("%5.20lf ", m(i,j));
        o << "\n";
    }
	
    return o;
}

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

Matrix& Matrix::operator / (Matrix &m) {
	if (this->n_column != m.n_row) {
		cout << "Matrix div: error in n_column/n_row\n";
		exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = &m.inv();
	
	return (*this) * (*m_aux);
}

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

Matrix& Matrix::transpose() {
	Matrix *m_aux = new Matrix(this->n_column, this->n_row);
	
	for(int i = 1; i <= this->n_row; i++) {
		for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(j,i) = (*this)(i,j);
		}
	}
	
	return *m_aux;
}

Matrix& zeros(const int n_row, const int n_column) {
	Matrix *m_aux = new Matrix(n_row, n_column);
	
	for(int i = 1; i <= n_row; i++) {
		for(int j = 1; j <= n_column; j++) {
			(*m_aux)(i,j) = 0;
		}
	}
	
	return (*m_aux);
}

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