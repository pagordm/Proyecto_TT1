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
#include "..\include\global.hpp"
#include "..\include\accelharmonic.hpp"
#include "..\include\eqnequinox.hpp"
#include "..\include\LTC.hpp"
#include "..\include\JPL_Eph_DE430.hpp"
#include "..\include\nutmatrix.hpp"
#include "..\include\PoleMatrix.hpp"
#include "..\include\PrecMatrix.hpp"
#include "..\include\gmst.hpp"
#include "..\include\gast.hpp"
#include "..\include\measupdate.hpp"
#include "..\include\g_accelharmonic.hpp"
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

	Matrix B = eye(f);
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
	Matrix C = eye(f);
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
	double result = Frac(x);
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
	Matrix s = eye(15);
	std::tuple<double, double, double, double, double, double, double, double, double> result = iers(eopdata, 37670);
	

	_assert(fabs(std::get<0>(result) - (-1.338037278494208e-07)) < 1e-10);
	_assert(fabs(std::get<1>(result) - 1.058353113998928e-06) < 1e-10);
	_assert(fabs(std::get<2>(result) - 0.030535300000000) < 1e-10);

	return 0;
}

int legendre_01() {
	auto [P, dP] = Legendre(1, 2, 1.0);
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
	Matrix r=eye(3);
	Matrix v(3);
	v(1)=1; v(2)=2; v(3)=3;
	Matrix result = TimeUpdate(r, v, 1.0);
	Matrix expected(1);
	expected(1)=15;

	_assert(m_equals(result, expected.transpose(), 1e-10));
	return 0;
}

int harmonic_01() {
	double a = -4.783104375562398e+09;
	Matrix r(3,3);
	r(1,1)=1; r(1,2)=2; r(1,3)=3;
	r(2,1)=4; r(2,2)=5; r(2,3)=6;
	r(3,1)=7; r(3,2)=8; r(3,3)=9;
	Matrix E(3);
	E(1)=10; E(2)=10; E(3)=10;
	Matrix expected(3);
	expected(1)=-3.851815228199548e+10; expected(2)=-4.592548925930230e+10; expected(3)=-5.333282623660913e+10;
	Matrix result = AccelHarmonic(E.transpose(),r, 1, 10);
	_assert(m_equals(expected.transpose(), result, 1e-10));
	return 0;
}
int eqnequinox_01() {

	double expected = 2.923587975442177e-05; //From matlab
	double result = EqnEquinox(10.0);

	//cout << result << ", expected=" << expected <<endl;

	_assert(fabs(expected-result) < 1e-9);

	return 0;
}

int jpl_01() {
	Matrix a(3);
	a(1)=8.381066621167137e+10; a(2)=-6.530025356227299e+10; a(3)=-2.340059142315123e+10;
	Matrix b(3);
	b(1)=-1.525487109073703e+10; b(2)=-1.101196198394985e+11; b(3)=-4.101455778214513e+10;
	Matrix c(3);
	c(1)=-9.244637046666153e+10; c(2)=1.064127619878867e+11; c(3)=4.613787655122181e+10;

	auto [ar, br, cr, d, e, f, g, h, i, j, k] = JPL_Eph_DE430(49746.1);

	_assert(m_equals(a.transpose(), ar, 1e-4));
	_assert(m_equals(b.transpose(), br, 1e-4));
	_assert(m_equals(c.transpose(), cr, 1e-4));

	return 0;
}

int ltc_01() {

	Matrix r(3,3);
	r(1,1)=0.157745694143248; r(1,2)=-0.987479769908865; r(1,3)=0;
	r(2,1)=-0.975661240821901; r(2,2)=-0.155857734378009; r(2,3)=0.154251449887584;
	r(3,1)=-0.152320186243100; r(3,2)=-0.024332502035119; r(3,3)=-0.988031624092862;

	Matrix result = LTC(3.3, 30.0);

	_assert(m_equals(r, result, 1e-10));
	return 0;

}

int nutmatrix_01() {

	double Mjd_TT = 54321.0;
	Matrix r(3,3);
	r(1,1)=0.999999999333437; r(1,2)=-0.000033499358637; r(1,3)=-0.000014523062858;
	r(2,1)=0.000033498769882; r(2,2)=0.999999998617265; r(2,3)=-0.000040537679027;
	r(3,1)=0.000014524420824; r(3,2)=0.000040537192496; r(3,3)=0.999999999072889;
	Matrix result = NutMatrix(Mjd_TT);

	_assert(m_equals(result, r, 1e-9));
	return 0;
}

int polemat_01() {

	Matrix r(3,3);
	r(1,1)=-0.839071529076452; r(1,2)=0.295958969093304; r(1,3)=0.456472625363814;
	r(2,1)=0; r(2,2)=-0.839071529076452; r(2,3)=0.544021110889370;
	r(3,1)=0.544021110889370; r(3,2)=0.456472625363814; r(3,3)=0.704041030906696;

	Matrix result = PoleMatrix(10.0, 10.0);

	_assert(m_equals(r, result, 1e-10));
	return 0;
}

int precmat_01() {
	Matrix r(3,3);
	r(1,1)=0.783082323064378; r(1,2)=0.571469950029840; r(1,3)=0.245365383697434;
	r(2,1)=-0.571567748001079; r(2,2)=0.816815000669506; r(2,3)=-0.078253205213914;
	r(3,1)=-0.245137481322363; r(3,2)=-0.078964238071218; r(3,3)=0.966267180626953;
	Matrix result = PrecMatrix(1e6, 50);
	
	_assert(m_equals(r, result, 1e-10));
	return 0;
}

int gmst_01() {
	double expected = 1.145236060990417;
	double result = gmst(10.0);

	_assert(fabs(expected-result) < 1e-10);
	return 0;
}

int gast_01() {
	double expected=1.145265296870172;

	double result = gast(10.0);
	//cout << fabs(expected-result) << endl;
	_assert(fabs(expected-result) < 1e-9);
	return 0;
}

int measupdate_01() {
	Matrix Y(6,1);
	Y(1,1)=7101576.98989518;
	Y(2,1)=1295199.87126989;
	Y(3,1)=12739.282331058;
	Y(4,1)=576.004647735755;
	Y(5,1)=-3084.62203921229;
	Y(6,1)=-6736.0259467681;
	double obs = 3.196905628244;
	double Azim = 3.19766548246716;
	double sigma_az = 0.00039095375244673;
	Matrix dAdY(6);
	dAdY(1)=6.8011430799027e-08; dAdY(2)=-3.73341445315677e-07; dAdY(3)=1.98045516802789e-08; dAdY(4)=0; dAdY(5)=0; dAdY(6)=0;
	Matrix P(6,6);
	P(1,1)=16662.1344187737; P(1,2)=-5939.90058143741;P(1,3)=8994.78627408343; P(1,4)=50.9584525113159;   P(1,5)=-14.0175161399365;  P(1,6)=23.634816671666;
	P(2,1)=-5939.90058143741;P(2,2)=25083.2687156787; P(2,3)=-1608.29906724942;P(2,4)=-5.16755530262711;  P(2,5)=61.0478532927088;   P(2,6)=-27.7748186275101;
	P(3,1)=8994.78627408342; P(3,2)=-1608.29906724941;P(3,3)=6462.95085906703; P(3,4)=28.1143150411363;   P(3,5)=-4.40165366994466;  P(3,6)=16.752629123315;
	P(4,1)=50.9584525113159; P(4,2)=-5.16755530262713;P(4,3)=28.1143150411364; P(4,4)=0.185799761897259;  P(4,5)=-0.0201977467654655;P(4,6)=0.0703233725080075;
	P(5,1)=-14.0175161399365;P(5,2)=61.0478532927088; P(5,3)=-4.4016536699447; P(5,4)=-0.0201977467654654;P(5,5)=0.156251627321938;  P(5,6)=-0.0770413411227855;
	P(6,1)=23.6348166716659; P(6,2)=-27.77481862751;  P(6,3)=16.752629123315;  P(6,4)=0.0703233725080072; P(6,5)=-0.0770413411227851;P(6,6)=0.086754864868731;

	Matrix expectedP(6,6);
	expectedP(1,1)=16582.6959715575;	expectedP(1,2)=-5719.28825199721;	expectedP(1,3)=8964.61806949935;	expectedP(1,4)=50.8244747796258;	expectedP(1,5)=-13.4810430619866;	expectedP(1,6)=23.3577425895417;	
	expectedP(2,1)=-5719.28825199721;	expectedP(2,2)=24470.5956133384;	expectedP(2,3)=-1524.51749646659;	expectedP(2,4)=-4.79547930520971;	expectedP(2,5)=59.5579881259538;	expectedP(2,6)=-27.0053428783355;	
	expectedP(3,1)=8964.61806949935;	expectedP(3,2)=-1524.51749646657;	expectedP(3,3)=6451.49393109675;	expectedP(3,4)=28.0634345448404;	expectedP(3,5)=-4.1979181975556;	expectedP(3,6)=16.6474051683235;	
	expectedP(4,1)=50.8244747796259;	expectedP(4,2)=-4.79547930520973;	expectedP(4,3)=28.0634345448405;	expectedP(4,4)=0.185573800373364;	expectedP(4,5)=-0.0192929525680047;	expectedP(4,6)=0.0698560703598162;	
	expectedP(5,1)=-13.4810430619866;	expectedP(5,2)=59.5579881259538;	expectedP(5,3)=-4.19791819755564;	expectedP(5,4)=-0.0192929525680046;	expectedP(5,5)=0.15262865414647;	expectedP(5,6)=-0.0751701717977754;	
	expectedP(6,1)=23.3577425895416;	expectedP(6,2)=-27.0053428783354;	expectedP(6,3)=16.6474051683235;	expectedP(6,4)=0.0698560703598159;	expectedP(6,5)=-0.075170171797775;	expectedP(6,6)=0.0857884556591442;	
	Matrix expectedY(6,1);
	expectedY(1,1)=7101559.88526246;	
	expectedY(2,1)=1295247.37336744;	
	expectedY(3,1)=12732.7865336534;	
	expectedY(4,1)=575.975799741193;	
	expectedY(5,1)=-3084.50652619203;	
	expectedY(6,1)=-6736.08560617204;	
	Matrix K(6,1);
	K(1,1)=22510.4134439327;
	K(2,1)=-62514.7509871813;
	K(3,1)=8548.74159612631;
	K(4,1)=37.9651697422456;
	K(5,1)=-152.019975331698;
	K(6,1)=78.5142756660613;
	auto [newK, newY, newP] = MeasUpdate(Y, obs, Azim, sigma_az, dAdY, P, 6.0);
	_assert(m_equals(K, newK, 1e-7));
	_assert(m_equals(expectedY, newY, 1e-7));
	_assert(m_equals(expectedP, newP, 1e-7));
	return 0;
}

int g_accelharmonic_01() {
	Matrix r(3,1);
	r(1,1)=7101800.90695315;
	r(2,1)=1293997.58115302;
	r(3,1)=10114.014948955;
	Matrix U(3,3);
	U(1,1)=-0.984320311904791; U(1,2)=0.17638970840918; U(1,3)=-0.000440838949610109;
	U(2,1)=-0.176389673507182; U(2,2)=-0.984320409906027; U(2,3)=-0.000117142904888635;
	U(3,1)=-0.000454589578418276; U(3,2)=-3.75467022865179e-05; U(3,3)=0.999999895969275;
	int n_max = 20;
	int m_max = 20;
	Matrix expected(3,3);
	expected(1,1)=2.02233500257165e-06; expected(1,2)=5.61803303433805e-07; expected(1,3)=4.39856240319614e-09;
	expected(2,1)=5.61803301435404e-07; expected(2,2)=-9.58631634517815e-07; expected(2,3)=8.05634892131479e-10;
	expected(3,1)=4.39855909334375e-09; expected(3,2)=8.0563404905587e-10; expected(3,3)=-1.06370336962723e-06;

	Matrix result = G_AccelHarmonic(r, U, n_max, m_max);

	_assert(m_equals(expected, result, 1e-7));
	return 0;

}

int all_tests()
{
	eop19620101(10); // c = 21413
	GGM03S();
	DE430Coeff();

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
	_verify(harmonic_01);
	_verify(jpl_01);
	_verify(eqnequinox_01);
	_verify(nutmatrix_01);
	_verify(polemat_01);
	_verify(precmat_01);
	_verify(gmst_01);
	_verify(gast_01);
	_verify(measupdate_01);
	_verify(g_accelharmonic_01);
	
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
