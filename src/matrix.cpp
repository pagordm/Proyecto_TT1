#include "..\include\matrix.h"

/**
 * @brief Construct a new Matrix:: Matrix object
 * 
 */
Matrix::Matrix() {
	this->n_column=0;
	this->n_row=0;
	this->data=nullptr;
}

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
 * @brief Construct a new Matrix:: vector object
 * 
 * @param v_size the size of the vector
 */
Matrix::Matrix(const int v_size) {
    if (v_size<=0) {
		cout << "Vector create: error in v_size\n";
        exit(EXIT_FAILURE);
	}
	
	this->n_row = 1;
	this->n_column = v_size;
	this->data = (double **) malloc(sizeof(double *));
	
    if (this->data == NULL) {
		cout << "Matrix create: error in data\n";
        exit(EXIT_FAILURE);
	}
	
	this->data[0] = (double *) calloc(v_size, sizeof(double));
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
		cout << "Matrix get: error in row/column: " << row << "," << column << endl;
        exit(EXIT_FAILURE);
	}
	
	return this->data[row - 1][column - 1];
}

/**
 * @brief Accesses a element of the vector.
 * 
 * @param row Row of the element
 * @param column Column of the element
 * @return double& A reference to the element in the row and columns
 */
double& Matrix::operator () (const int n) {
	if (n<=0 || n>this->n_row*this->n_column) {
		cout << "Matrix get: error in row/column: " << n << endl;
        exit(EXIT_FAILURE);
	}
	
	return this->data[(n - 1)/this->n_column][(n - 1)%this->n_column];
}

/**
 * @brief Assignment operator. Copies all values of matrix m into this. Requires both matrices to have the same shape.
 * 
 * @param m The matrix to copy data from
 * @return Matrix& A reference to this, to be used in chaining assignments.
 */
Matrix& Matrix::operator = (Matrix &m) {
	//Recreate the matrix with the same shape as m
	this->n_row = m.n_row;
	this->n_column = m.n_column;

	this->data = (double **) malloc(m.n_row*sizeof(double *));
	
    if (this->data == NULL) {
		cout << "Matrix assignment: error in data\n";
        exit(EXIT_FAILURE);
	}
	
	for(int i = 0; i < m.n_row; i++) {
		this->data[i] = (double *) malloc(m.n_column*sizeof(double));
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
		cout << "Matrix mul: error in n_column/n_row" << this->n_column << " " << m.n_row << endl;
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
	
	Matrix *m_aux = &inv(m);
	
	return (*this) * (*m_aux);
}
/**
 * @brief Calculates the inverse of this matrix
 * 
 * @return Matrix& Reference to the inverted matrix
 */
Matrix& inv(Matrix &m) {
	if (m.n_row != m.n_column) {
		cout << "Matrix inv: error in n_row/n_column\n";
		exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = &eye(m.n_row, m.n_column);
	Matrix m_copy = m;
	for(int i = 1; i <= m_copy.n_row; i++) {
		double pivot = m_copy(i,i);

		if (pivot==0) {
			cout << "Matrix inv: singular matrix not invertible\n";
			exit(EXIT_FAILURE);
		}
		
		for(int j = 1; j <= m_copy.n_column; j++) {
			m_copy(i,j) /= pivot;
			(*m_aux)(i,j) /= pivot;
		}
		
		for(int k = 1; k <= m_copy.n_row; k++) {
			if (k != i) {
				double factor = m_copy(k,i);
				
				for(int j = 1; j <= m_copy.n_column; j++) {
					m_copy(k,j) -= factor * m_copy(i,j);
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
 * @brief Generates a vector with all zeros.
 * 
 * @param n number of elements of the vector
 * @return Matrix& A reference to the vector
 */
Matrix& zeros(const int n) {
	Matrix *m_aux = new Matrix(n);
	
	for(int i = 1; i <= n; i++) {
		(*m_aux)(i) = 0;
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

/**
 * @brief Subtracts a matrix element-wise by a double
 * 
 * @param n Double to subtract the matrix by.
 * @return Matrix& Reference to the resulting matrix.
 */
Matrix& Matrix::operator - (const double n) {
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);

	for(int i = 1; i <= this->n_row; i++) {
		for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j)-n;
		}
	}
	return *m_aux;

}
/**
 * @brief Calculates the dot product between two vectors
 * 
 * @param m1 first vector of the dot product
 * @param m2 second vector of the dot product
 * @return double reference to the dot product
 */
double dot(Matrix &m1, Matrix &m2) {
	if (m1.n_row!=1 || m2.n_row!=1 || m1.n_column!=m2.n_column) {
		cout << "Vector dot: error in arguments" << endl;
		exit(EXIT_FAILURE);
	}
	double sum = 0;
	for(int i = 1; i <=m1.n_column; i++) {
		sum+=m1(i)*m2(i);
	}

	return sum;
}
/**
 * @brief Calculates the norm of a vector
 * 
 * @param m1 vector to calculate the norm of
 * @return double the norm
 */
double norm(Matrix &m1) {
	double sum = 0.0;
	for(int i = 1; i <=m1.n_column; i++) {
		sum+=m1(i)*m1(i);
	}
	return sqrt(sum);
}

/**
 * @brief Calculates the cross product between two vectors
 * 
 * @param m1 first vector of the cross product
 * @param m2 second vector of the cross product
 * @return Matrix& reference to the cross product matrix
 */
Matrix& cross(Matrix &m1, Matrix &m2) {
	if (m1.n_column!=3 || m2.n_column!=3 || m1.n_row!=1 || m2.n_row!=1) {
		cout << "Vector cross: error in arguments" << endl;
		exit(EXIT_FAILURE);
	}
	Matrix *m_aux = new Matrix(3);

	(*m_aux)(1) = m1(2)*m2(3) - m1(3)*m2(2);
	(*m_aux)(2) = m1(3)*m2(1) - m1(1)*m2(3);
	(*m_aux)(3) = m1(1)*m2(2) - m1(2)*m2(1);
	
	return *m_aux;
	
}

/**
 * @brief Extracts a subvector from a vector between specified indices
 * 
 * @param start starting index of extraction (inclusive)
 * @param end ending index of extraction (inclusive) 
 * @return Matrix& reference to the extracted vector
 */
Matrix& Matrix::extract_vector(const int start, const int end) {
    if (this->n_row != 1 || start < 1 || end > this->n_column || start > end) {
        cout << "Extract vector: error in arguments" << endl;
        exit(EXIT_FAILURE);
    }

    int new_size = end - start + 1;
    Matrix *result = new Matrix(new_size);

    for (int i = 0; i < new_size; i++) {
        (*result)(i+1) = (*this)(start+i);
    }

    return *result;
}

/**
 * @brief Extracts a row from a matrix
 * 
 * @param n the row to extract
 * @return Matrix& reference to the extracted row as a vector
 */
Matrix& Matrix::extract_row(const int n) {
    if (n < 1 || n > this->n_row) {
        cout << "Extract row: error in arguments" << endl;
        exit(EXIT_FAILURE);
    }

    Matrix *result = new Matrix(n_column);
    
    for (int j = 0; j < n_column; j++) {
        (*result)(j+1) = (*this)(n, j+1);
    }

    return *result;
}

/**
 * @brief Extracts a column from a matrix
 * 
 * @param n the column to extract
 * @return Matrix& reference to the extracted column as a vector
 */
Matrix& Matrix::extract_column(const int n) {
    if (n < 1 || n > this->n_column) {
        cout << "Extract column: error in arguments" << endl;
        exit(EXIT_FAILURE);
    }
	Matrix *result = new Matrix(this->n_row);

	for (int i = 0; i < this->n_row; i++) {
		(*result)(i+1) = (*this)(i+1, n);
	}

	return *result;
}

/**
 * @brief Concatenates two vectors into a single vector
 * 
 * @param m1 first vector to concatenate
 * @param m2 second vector to concatenate
 * @return Matrix& reference to the concatenated vector
 */
Matrix& union_vector(Matrix &m1, Matrix &m2) {
    if (m1.n_row != 1 || m2.n_row != 1 || m1.n_column != m2.n_column) {
        cout << "Union vector: error in arguments" << endl;
        exit(EXIT_FAILURE);
    }

    Matrix *result = new Matrix(m1.n_column + m2.n_column);

    for (int i = 0; i < m1.n_column; i++) {
        (*result)(i+1) = m1(i+1);
	}

    for (int i = 0; i < m2.n_column; i++) {
        (*result)(i+m1.n_column+1) = m2(i+1);
    }

    return *result;
}

/**
 * @brief Assigns a row to a matrix
 * 
 * @param n the row to assign
 * @param m the matrix to assign
 * @return Matrix& reference to the matrix
 */
Matrix& Matrix::assign_row(const int n, Matrix &m) {
	if (n < 1 || n > this->n_row || m.n_row != 1 || m.n_column != this->n_column) {
		cout << "Assign row: error in arguments" << endl;
		exit(EXIT_FAILURE);
	}
	
	for (int i = 0; i < this->n_column; i++) {
		(*this)(n, i+1) = m(i+1);
	}

	return *this;
}

/**
 * @brief Assigns a column to a matrix
 * 
 * @param n the column to assign
 * @param m the matrix to assign (must be a vector)
 * @return Matrix& reference to the matrix
 */
Matrix& Matrix::assign_column(const int n, Matrix &m) {
	if (n < 1 || n > this->n_column || m.n_row != 1 || m.n_column != this->n_row) {
		cout << "Assign column: error in arguments" << endl;
		exit(EXIT_FAILURE);
	}
	
	for (int i = 0; i < this->n_column; i++) {
		(*this)(i+1, n) = m(i+1);
	}

	return *this;
}


