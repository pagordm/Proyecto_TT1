#define _USE_MATH_DEFINES
#include "..\include\matrix.h"
#include "..\include\R_x.hpp"
#include "..\include\R_y.hpp"
#include "..\include\R_z.hpp"
#include "..\include\position.hpp"
#include "..\include\accelpointmass.hpp"
#include "..\include\cheb3d.hpp"
#include "..\include\sign_.hpp"
#include "..\include\eccanom.hpp"
#include "..\include\frac.hpp"
#include "..\include\meanobliquity.hpp"
#include "..\include\mjday.hpp"
#include "..\include\mjday_tdb.hpp"
#include "..\include\timediff.hpp"
#include "..\include\nutangles.hpp"
#include "..\include\timeupdate.hpp"
#include "..\include\azelpa.hpp"
#include "..\include\iers.hpp"
#include "..\include\legendre.hpp"
#include <tuple>
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

int R_x_01() {
	Matrix R = R_x(M_PI/2);
	Matrix expected(3,3);
	expected(1,1)=1; expected(1,2)=0; expected(1,3)=0;
	expected(2,1)=0; expected(2,2)=0; expected(2,3)=1;
	expected(3,1)=0; expected(3,2)=-1; expected(3,3)=0;

	_assert(m_equals(R, expected, 1e-10));
	return 0;
}

int R_y_01() {
	Matrix R = R_y(M_PI/2);
	Matrix expected(3,3);
	expected(1,1)=0; expected(1,2)=0; expected(1,3)=-1;
	expected(2,1)=0; expected(2,2)=1; expected(2,3)=0;
	expected(3,1)=1; expected(3,2)=0; expected(3,3)=0;

	_assert(m_equals(R, expected, 1e-10));
	return 0;
}

int R_z_01() {
	Matrix R = R_z(M_PI/2);
	Matrix expected(3,3);
	expected(1,1)=0; expected(1,2)=1; expected(1,3)=0;
	expected(2,1)=-1; expected(2,2)=0; expected(2,3)=0;
	expected(3,1)=0; expected(3,2)=0; expected(3,3)=1;

	_assert(m_equals(R, expected, 1e-10));
	return 0;
}

int position_01() {
	Matrix p = position(0, 0, 0);
	Matrix expected(3);
	expected(1)=6378136.3; expected(2)=0; expected(3)=0;

	_assert(m_equals(p, expected, 1e-10));
	return 0;
}

int accelpointmass_01() {
	Matrix r(3);
	r(1)=1; r(2)=2; r(3)=3;
	Matrix s(3);
	s(1)=4; s(2)=5; s(3)=6;
	Matrix a = accelpointmass(r, s, 1);
	Matrix expected(3);
	expected(1)=0.015463313357364; expected(2)=0.013983305870875; expected(3)=0.012503298384387; //results from MATLAB

	_assert(m_equals(a, expected, 1e-10));
	return 0;
}

int cheb3d_01() {
	Matrix Cx(2);
	Cx(1)=1; Cx(2)=2;
	Matrix Cy(2);
	Cy(1)=2; Cy(2)=3;
	Matrix Cz(2);
	Cz(1)=3; Cz(2)=4;
	
	Matrix result = cheb3d(1, 2, 0, 2, Cx, Cy, Cz);
	Matrix expected(3);
	expected(1)=1; expected(2)=2; expected(3)=3;

	_assert(m_equals(result, expected, 1e-10));
	return 0;
}

int sign_01() {
	double a = -5.0;
	double b = -10.0;
	double result = sign_(a, b);
	double expected = -5.0;

	_assert(fabs(result - expected) < 1e-10);
	return 0;
}

int eccanom_01() {
	double M = 0.5;
	double e = 0.0;
	double result = EccAnom(M, e);
	double expected = 0.5;

	_assert(fabs(result - expected) < 1e-10);
	return 0;
}

int frac_01() {
	double x = 1.5;
	double result = frac(x);
	double expected = 0.5;

	_assert(fabs(result - expected) < 1e-10);
	return 0;
}

int meanobliquity_01() {
	double Mjd_TT = 54321.0;
	double result = MeanObliquity(Mjd_TT);
	double expected = 0.409075551101389; //results from MATLAB

	_assert(fabs(result - expected) < 1e-10);
	return 0;
}

int mjday_01() {
	int year = 2025;
	int month = 1;
	int day = 1;
	int hour = 1;
	int minute = 0;
	double second = 0;
	double result = Mjday(year, month, day, hour, minute, second);
	double expected = 6.067604166666651e04; //results from MATLAB

	_assert(fabs(result - expected) < 1e-10);
	return 0;
}

int mjday_tdb_01() {
	double Mjd_TT = 6.067604166666651e04;
	double result = Mjday_TDB(Mjd_TT);
	double expected = 6.067604166666546e04; //results from MATLAB

	_assert(fabs(result - expected) < 1e-10);
	return 0;
}

int timediff_01() {
	double Mjd_TT = 6.067604166666651e04;
	double Mjd_UTC = 7.0e04;
	std::tuple<double, double, double, double, double> result = timediff(Mjd_TT, Mjd_UTC);
	std::tuple<double, double, double, double, double> expected = std::make_tuple(-9.323958333333489e+03, -69981, -9.304958333333489e+03, 7.003218399999999e+04, 69981);
	_assert(fabs(std::get<0>(result) - std::get<0>(expected)) < 1e-10);
	_assert(fabs(std::get<1>(result) - std::get<1>(expected)) < 1e-10);
	_assert(fabs(std::get<2>(result) - std::get<2>(expected)) < 1e-10);
	_assert(fabs(std::get<3>(result) - std::get<3>(expected)) < 1e-10);
	_assert(fabs(std::get<4>(result) - std::get<4>(expected)) < 1e-10);
	return 0;
}

int azelpa_01() {
	Matrix r(3);
	r(1)=1; r(2)=2; r(3)=3;
	double expected = 0.463647609000806; //results from MATLAB
	double expected2 = 0.930274014115472;
	std::tuple<double, double, Matrix, Matrix> result = AzElPa(r);

	_assert(fabs(std::get<0>(result) - expected) < 1e-10);
	_assert(fabs(std::get<1>(result) - expected2) < 1e-10);
	return 0;
}

int iers_01() {
	Matrix s = eye(15, 15);
	std::tuple<double, double, double, double, double, double, double, double, double> result = iers(s, 0.0);
	

	_assert(fabs(std::get<0>(result) - 0.0) < 1e-10);
	_assert(fabs(std::get<1>(result) - 0.0) < 1e-10);
	_assert(fabs(std::get<2>(result) - 0.0) < 1e-10);

	return 0;
}

int legendre_01() {
	std::tuple<Matrix, Matrix> result = legendre(1, 2, 1.0);
	Matrix P = std::get<0>(result);
	Matrix dP = std::get<1>(result);
	Matrix expected(2,3);
	expected(1,1) = 1.0;    expected(1,2) = 0.0;       expected(1,3) = 0.0;
	expected(2,1) = 1.457470498782296; expected(2,2) = 0.935831045210238;    expected(2,3) = 0.0;
	_assert(m_equals(P, expected, 1e-10));
	return 0;
}

int nutangles_01() {
	double Mjd_TT = 6.067604166666651e04;
	double expected = 9.723503682287755e-07;
	double expected2 = 4.120518762261807e-05;
	std::tuple<double, double> result = NutAngles(Mjd_TT);

	_assert(fabs(std::get<0>(result) - expected) < 1e-10);
	_assert(fabs(std::get<1>(result) - expected2) < 1e-10);
	return 0;
}

int timeupdate_01() {
	Matrix r=eye(3,3);
	Matrix v(3);
	v(1)=1; v(2)=2; v(3)=3;
	Matrix result = TimeUpdate(r, v, 1.0);
	Matrix expected(1);
	expected(1)=15;

	_assert(m_equals(result, expected.transpose(), 1e-10));
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
	_verify(R_x_01);
	_verify(R_y_01);
	_verify(R_z_01);
	_verify(position_01);
	_verify(accelpointmass_01);
	_verify(cheb3d_01);
	_verify(sign_01);
	_verify(eccanom_01);
	_verify(frac_01);
	_verify(meanobliquity_01);
	_verify(mjday_01);
	_verify(mjday_tdb_01);
	_verify(timediff_01);
	_verify(azelpa_01);
	_verify(iers_01);
	_verify(legendre_01);
	_verify(nutangles_01);
	_verify(timeupdate_01);

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
