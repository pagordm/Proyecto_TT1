#include "..\include\matrix.h"
#include <cstdio>
#include <cmath>

int tests_run = 0;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

int m_equals(Matrix A, Matrix B, double p) {
	if (A.n_row != B.n_row || A.n_column != B.n_column)
		return 0;
	else
		for(int i = 1; i <= A.n_row; i++)
			for(int j = 1; j <= A.n_column; j++)
				if(fabs(A(i,j)-B(i,j)) > p) {
					printf("%2.20lf %2.20lf\n",A(i,j),B(i,j));
					return 0;
				}
	
	return 1;
}

int m_sum_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) =  0; B(1,3) = 0; B(1,4) = 0;
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; B(3,4) = 2;
	
	Matrix C(f, c);
	C(1,1) = 2; C(1,2) =  2; C(1,3) = 8; C(1,4) = 0;
	C(2,1) = 8; C(2,2) = -3; C(2,3) = 1; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = -2; C(3,3) = 0; C(3,4) = 7;
	
	Matrix R = A + B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_sub_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) = 0; B(1,3) = 0; B(1,4) = 0;
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; B(3,4) = 2;
	
	Matrix C(f, c);
	C(1,1) = -2; C(1,2) = 2; C(1,3) = 8; C(1,4) = 0;
	C(2,1) = -6; C(2,2) = 1; C(2,3) = -1; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = 4; C(3,3) = 0; C(3,4) = 3;
	
	Matrix R = A - B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_zeros_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 0; A(1,3) = 0; A(1,4) = 0;
	A(2,1) = 0; A(2,2) = 0; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 0; A(3,4) = 0;
	
	Matrix B = zeros(3, 4);
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}

int m_inv_01() {
	int f = 2;
	int c = 2;
	Matrix A(f, c);
	A(1,1) = 1; A(1,2) = 2;
	A(2,1) = 2; A(2,2) = 3;

	Matrix B(f, c); //Matriz inversa de A
	B(1,1) = -3; B(1,2) = 2;
	B(2,1) = 2; B(2,2) = -1;
	_assert(m_equals(inv(A), B, 1e-10));
	return 0;
}

int m_eye_01() {
	int f = 3;
	int c = 3;
	Matrix A(f, c);
	A(1,1) = 1; A(1,2) = 0; A(1,3) = 0;
	A(2,1) = 0; A(2,2) = 1; A(2,3) = 0;
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 1;

	Matrix B = eye(f, c);
	_assert(m_equals(A, B, 1e-10));

	return 0;
}

int m_mul_01() {
	int f = 3;
	int c = 3;

	Matrix A(f, c);
	A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
	A(2,1) = 4; A(2,2) = 5; A(2,3) = 6;
	A(3,1) = 7; A(3,2) = 8; A(3,3) = 9;

	Matrix B(f, c);
	B(1,1) = 9; B(1,2) = 8; B(1,3) = 7;
	B(2,1) = 6; B(2,2) = 5; B(2,3) = 4;
	B(3,1) = 3; B(3,2) = 2; B(3,3) = 1;

	Matrix C(f, c);
	C(1,1) = 30; C(1,2) = 24; C(1,3) = 18;
	C(2,1) = 84; C(2,2) = 69; C(2,3) = 54;
	C(3,1) = 138; C(3,2) = 114; C(3,3) = 90;

	Matrix R = A * B;
	_assert(m_equals(C, R, 1e-10));

	return 0;
}
int m_div_01() {
	int f = 2;
	int c = 2;
	Matrix A(f, c);
	A(1,1) = 1; A(1,2) = 2;
	A(2,1) = 2; A(2,2) = 3;
	Matrix B(f, c);
	B(1,1) = 1; B(1,2) = 2;
	B(2,1) = 2; B(2,2) = 3;
	Matrix C = eye(f, c);
	_assert(m_equals(A / B, C, 1e-10));
	return 0;
}

int m_div_02() {
	int f = 2;
	int c = 2;
	Matrix A(f,c);
	A(1,1) = 2; A(1,2) = 4;
	A(2,1) = 4; A(2,2) = 6;
	Matrix B(f, c);
	B(1,1) = 1; B(1,2) = 2;
	B(2,1) = 2; B(2,2) = 3;
	Matrix C = A/2.0;

	_assert(m_equals(B, C, 1e-10));
	return 0;
}

int m_eq_01() {
	int f = 2;
	int c = 2;
	Matrix A(f,c);
	A(1,1) = 1; A(1,2) = 2;
	A(2,1) = 2; A(2,2) = 3;
	Matrix B = A;

	_assert(m_equals(A, B, 1e-10));
	return 0;
}
int m_eq_02() {
	int f = 2;
	int c = 2;
	Matrix A(f,c);
	A(1,1) = 2; A(1,2) = 4;
	A(2,1) = 4; A(2,2) = 6;
	Matrix C(f, c);
	C(1,1) = 1; C(1,2) = 2;
	C(2,1) = 2; C(2,2) = 3;
	Matrix B(f,c);
	B = A/2.0;

	_assert(m_equals(B, C, 1e-10));
	_assert(m_equals(A/2.0, B, 1e-10));
	return 0;
}

int m_mul_02() {
	int f = 2;
	int c = 2;
	Matrix A(f,c);
	A(1,1) = 2; A(1,2) = 4;
	A(2,1) = 4; A(2,2) = 6;
	Matrix C(f, c);
	C(1,1) = 1; C(1,2) = 2;
	C(2,1) = 2; C(2,2) = 3;

	_assert(m_equals(C*2.0, A, 1e-10));
	return 0;
}

int m_add_02() {
	int f = 2;
	int c = 2;
	Matrix A(f,c);
	A(1,1) = 2; A(1,2) = 3;
	A(2,1) = 3; A(2,2) = 4;
	Matrix C(f, c);
	C(1,1) = 1; C(1,2) = 2;
	C(2,1) = 2; C(2,2) = 3;

	_assert(m_equals(C+1.0, A, 1e-10));
	return 0;
}

int m_sub_02() {
	int f = 2;
	int c = 2;
	Matrix A(f,c);
	A(1,1) = 2; A(1,2) = 3;
	A(2,1) = 3; A(2,2) = 4;
	Matrix C(f, c);
	C(1,1) = 1; C(1,2) = 2;
	C(2,1) = 2; C(2,2) = 3;

	_assert(m_equals(A-1.0, C, 1e-10));
	return 0;
}

int m_zeros_02() {
    int f = 1;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 0; A(1,3) = 0; A(1,4) = 0;
	
	Matrix B = zeros(4);
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}

int m_dot_01() {
	Matrix v1(3);
	v1(1)=1; v1(2)=2; v1(3)=3;
	Matrix v2(3);
	v2(1)=1; v2(2)=2; v2(3)=3;

	_assert(fabs(dot(v1, v2)-14.0)<1e-10);

	return 0;
}

int m_norm_01() {
	Matrix v1(2);
	v1(1)=3; v1(2)=4;
	_assert(fabs(norm(v1)-5.0)<1e-10);
	return 0;
}

int m_cross_01() {
	Matrix v1(3);
	v1(1)=1; v1(2)=0; v1(3)=0;
	Matrix v2(3);
	v2(1)=0; v2(2)=1; v2(3)=0;
	
	Matrix C(3);
	C(1)=0; C(2)=0; C(3)=1;
	
	Matrix R = cross(v1,v2);
	_assert(m_equals(C, R, 1e-10));
	return 0;
}

int m_extract_vector_01() {
	Matrix v(3);
	v(1)=1; v(2)=2; v(3)=3;
	Matrix subv = v.extract_vector(2,3);
	Matrix subv_expected(2);
	subv_expected(1)=2; subv_expected(2)=3;

	_assert(m_equals(subv, subv_expected, 1e-10));
	return 0;
}

int m_extract_row_01() {
	Matrix A(3,3);
	A(1,1)=1; A(1,2)=2; A(1,3)=3;
	A(2,1)=4; A(2,2)=5; A(2,3)=6;
	A(3,1)=7; A(3,2)=8; A(3,3)=9;
	
	Matrix row = A.extract_row(2);
	Matrix row_expected(3);
	row_expected(1)=4; row_expected(2)=5; row_expected(3)=6;

	_assert(m_equals(row, row_expected, 1e-10));
	return 0;
}

int m_union_vector_01() {
	Matrix v1(1,3);
	v1(1,1)=1; v1(1,2)=2; v1(1,3)=3;
	Matrix v2(1,3);
	v2(1,1)=4; v2(1,2)=5; v2(1,3)=6;
	Matrix v3(1,6);
	v3(1,1)=1; v3(1,2)=2; v3(1,3)=3; v3(1,4)=4; v3(1,5)=5; v3(1,6)=6;

	_assert(m_equals(union_vector(v1, v2), v3, 1e-10));
	return 0;
}

int m_extract_column_01() {
	Matrix A(3,3);
	A(1,1)=1; A(1,2)=2; A(1,3)=3;
	A(2,1)=4; A(2,2)=5; A(2,3)=6;
	A(3,1)=7; A(3,2)=8; A(3,3)=9;
	Matrix column = A.extract_column(2);
	Matrix column_expected(3);
	column_expected(1)=2; column_expected(2)=5; column_expected(3)=8;

	_assert(m_equals(column, column_expected, 1e-10));
	return 0;
}

int m_assign_row_01() {
	Matrix A(3,3);
	A(1,1)=1; A(1,2)=2; A(1,3)=3;
	A(2,1)=4; A(2,2)=5; A(2,3)=6;
	A(3,1)=7; A(3,2)=8; A(3,3)=9;
	Matrix row(3);
	row(1)=10; row(2)=11; row(3)=12;
	A.assign_row(2, row);
	Matrix A_expected(3,3);
	A_expected(1,1)=1; A_expected(1,2)=2; A_expected(1,3)=3;
	A_expected(2,1)=10; A_expected(2,2)=11; A_expected(2,3)=12;
	A_expected(3,1)=7; A_expected(3,2)=8; A_expected(3,3)=9;

	_assert(m_equals(A, A_expected, 1e-10));
	return 0;
}

int m_assign_column_01() {
	Matrix A(3,3);
	A(1,1)=1; A(1,2)=2; A(1,3)=3;
	A(2,1)=4; A(2,2)=5; A(2,3)=6;
	A(3,1)=7; A(3,2)=8; A(3,3)=9;
	Matrix column(3);
	column(1)=10; column(2)=11; column(3)=12;
	A.assign_column(2, column);
	Matrix A_expected(3,3);
	A_expected(1,1)=1; A_expected(1,2)=10; A_expected(1,3)=3;
	A_expected(2,1)=4; A_expected(2,2)=11; A_expected(2,3)=6;
	A_expected(3,1)=7; A_expected(3,2)=12; A_expected(3,3)=9;

	_assert(m_equals(A, A_expected, 1e-10));
	return 0;
}

int all_tests()
{
    _verify(m_sum_01);
    _verify(m_sub_01);
    _verify(m_zeros_01);
	_verify(m_inv_01);
	_verify(m_eye_01);
	_verify(m_mul_01);
	_verify(m_div_01);
	_verify(m_div_02);
	_verify(m_eq_01);
	_verify(m_eq_02);
	_verify(m_mul_02);
	_verify(m_add_02);
	_verify(m_sub_02);
	_verify(m_zeros_02);
	_verify(m_dot_01);
	_verify(m_norm_01);
	_verify(m_cross_01);
	_verify(m_extract_vector_01);
	_verify(m_extract_row_01);
	_verify(m_union_vector_01);
	_verify(m_extract_column_01);
	_verify(m_assign_row_01);
	_verify(m_assign_column_01);
    return 0;
}


int main()
{
    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}
