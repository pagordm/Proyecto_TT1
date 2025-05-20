/**
 * @file tests.cpp
 * @author Pablo Gordillo Minchinela
 * @brief This file contains the tests of the project.
 * @date 2025-05-20
 * 
 * 
 */
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
#include "..\include\GHAMatrix.hpp"
#include "..\include\accel.hpp"
#include "..\include\vareqn.hpp"
#include "..\include\geodetic.hpp"
#include "..\include\angl.hpp"
#include "..\include\elements.hpp"
#include "..\include\gibbs.hpp"
#include "..\include\hgibbs.hpp"
#include "..\include\unit.hpp"
#include "..\include\DEInteg.hpp"
#include "..\include\anglesg.hpp"
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
	r(1)=6.221397628578685e+06; r(2)=2.867713779657379e+06; r(3)=3.006155985099489e+06;
	Matrix s(3);
	s(1)=9.229825172847661e+10; s(2)=-1.053751960790543e+11; s(3)=-4.568636722635329e+10;
	Matrix a = AccelPointMass(r.transpose(), s.transpose(), 1.327124400419394e+20);
	Matrix expected(3);
	expected(1)=-1.8685505934417e-07; expected(2)=-2.00332995883182e-07; expected(3)=-1.59993120755489e-07; //results from MATLAB
	//cout << "result: \n" << a << endl;
	//cout << "expected: \n" << expected.transpose() << endl; 
	_assert(m_equals(a, expected.transpose(), 1e-10));
	return 0;
}

int cheb3d_01() {
	Matrix Cx(10);
	Cx(1)=0.040415016726888; Cx(2)=-0.000694670662256884; Cx(3)=0.00025866083199944; Cx(4)=7.50904395359041e-07; Cx(5)=-5.93855864146682e-07; Cx(6)=8.16521467824851e-07; Cx(7)=-1.91863736412801e-07; Cx(8)=-2.07578684050435e-08; Cx(9)=7.56679389322787e-09; Cx(10)=1.37625072631299e-10;
	Matrix Cy(10);
	Cy(1)=0.429712722039214; Cy(2)=0.000564692852447144; Cy(3)=3.26749688165755e-05; Cy(4)=-1.23295138693416e-05; Cy(5)=1.2265751042446e-06; Cy(6)=-2.62336840978799e-07;Cy(7)=-6.30254512604998e-08;Cy(8)=1.70983120176071e-08;Cy(9)=9.51078697383379e-10;Cy(10)=-5.50393584452035e-10;
	Matrix Cz(10);
	Cz(1)=2151.02929751302; Cz(2)=0.920502077686378; Cz(3)=-0.000224459082502766; Cz(4)=7.84520821145058e-07; Cz(5)=2.41635748679829e-07; Cz(6)=-7.67349910644716e-07; Cz(7)=1.8113846680676e-07; Cz(8)=1.91612698581954e-08; Cz(9)=-7.06484757375274e-09; Cz(10)=-1.23958515454317e-10;
	
	Matrix result = cheb3d(49746.1107720813, 10, 49744, 49752, Cx, Cy, Cz);
	Matrix expected(3);
	expected(1)=0.0406001295440521; expected(2)=0.429415265875236; expected(3)=2150.59466362242;

	_assert(m_equals(result, expected, 1e-10));
	return 0;
}

int cheb3d_02() { //Usamos r_Earth de JPL
	Matrix Cx(13);
	Cx(1) = -103506598.42109;
	Cx(2) = -14927175.9535129;
	Cx(3) = 512386.261705712;
	Cx(4) = 11989.7464659321;
	Cx(5) = -231.321791175806;
	Cx(6) = -3.01145680139332;
	Cx(7) = 0.0557108554642511;
	Cx(8) = 0.000160763608359754;
	Cx(9) = 4.5930459903819e-05;
	Cx(10) = 2.65245498743318e-06;
	Cx(11) = -2.69798189927076e-06;
	Cx(12) = 5.19937618744697e-07;
	Cx(13) = 5.4216592061808e-08;
	Matrix Cy(13);
	Cy(1)=96607067.3707606; 
	Cy(2)=-13337710.4674864; 
	Cy(3)=-474105.009264047; 
	Cy(4)=11293.5192596881; 
	Cy(5)=193.380285507684; 
	Cy(6)=-3.43928289232978; 
	Cy(7)=-0.0344227898113099; 
	Cy(8)=0.00107593031980365; 
	Cy(9)=4.63124437216978e-06; 
	Cy(10)=-1.2434711625373e-05; 
	Cy(11)=2.21298188012278e-06; 
	Cy(12)=1.39759592836155e-07; 
	Cy(13)=-1.68301628426181e-07;
	Matrix Cz(13);
	Cz(1) = 41887005.1061263; 
	Cz(2) = -5782554.14940129; 
	Cz(3) = -205553.6189751; 
	Cz(4) = 4896.48746963361; 
	Cz(5) = 83.8420287547238; 
	Cz(6) = -1.49161452214215;
	Cz(7) = -0.0148800126199173; 
	Cz(8) = 0.000502341000673383; 
	Cz(9) = -2.65374413641448e-06; 
	Cz(10) = -6.64088213358266e-06; 
	Cz(11) = 1.33131703155108e-06; 
	Cz(12) = 4.77273410670398e-08;
	Cz(13) = -9.29357771479359e-08;
	
	Matrix result = cheb3d(49746.1119928785, 13, 49744.0, 49744.0+16.0, Cx, Cy, Cz);
	Matrix expected(3);
	expected(1)=-9.246987516899471e+07; expected(2)=1.063908284562656e+08; expected(3)=4.612874685964671e+07;

	_assert(m_equals(result, expected, 1e-6)); //por qué 1e-6? No podemos hacerlo mejor?
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
	std::tuple<double, double, double, double, double, double, double, double, double> result = IERS(eopdata, 37670);
	

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
	// cout << "a: " << a << endl;
	// cout << "ar: " << ar << endl;
	_assert(m_equals(a.transpose(), ar, 1e-10));
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

	double Mjd_TT = 49746.1107720731;
	Matrix r(3,3);
	r(1,1)=0.999999998058926; r(1,2)=-5.71651736108225e-05; r(1,3)=-2.47848948831635e-05;
	r(2,1)=5.71660437946361e-05; r(2,2)=0.99999999774966; r(2,3)=3.51101546042787e-05;
	r(3,1)=2.47828877493056e-05; r(3,2)=-3.51115713905226e-05; r(3,3)=0.999999999076493;
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
	r(1,1)=0.999999279427994; r(1,2)=0.001101012138263; r(1,3)=0.000478451421900;
	r(2,1)=-0.001101012138264; r(2,2)=0.999999393885917; r(2,3)=-0.000000263388277;
	r(3,1)=-0.000478451421898; r(3,2)=-0.000000263392736; r(3,3)=0.999999885542077;
	Matrix result = PrecMatrix(5.154450000000000e+04, 4.974611085861109e+04);
	
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

int ghamatrix_01() {
	Matrix expected(3,3);
	expected(1,1)=0.412804512414729; expected(1,2)=0.910819649837463; expected(1,3)=0;
	expected(2,1)=-0.910819649837463; expected(2,2)=0.412804512414729; expected(2,3)=0;
	expected(3,1)=0; expected(3,2)=0; expected(3,3)=1.0;

	Matrix result = GHAMatrix(10.0);

	_assert(m_equals(expected, result, 1e-9));

	return 0;
}

int accel_01() {
	double x = -543.476874884521;
	Matrix Y(6, 1);
	Y(1,1)=5720694.2260585;
	Y(2,1)=2687728.41425142;
	Y(3,1)=3483000.08675422;
	Y(4,1)=4371.83136151615;
	Y(5,1)=-1905.47309296258;
	Y(6,1)=-5698.58341612187;
	Matrix expected(6,1);
	expected(1,1)=4371.83136151615;
	expected(2,1)=-1905.47309296258;
	expected(3,1)=-5698.58341612187;
	expected(4,1)=-6.0654420426171;
	expected(5,1)=-2.84977703178268;
	expected(6,1)=-3.70232534578346;
	Matrix result = Accel(x, Y);

	_assert(m_equals(result, expected, 1e-3));
	return 0;	
}

int vareqn_01() {
	double x = 5.38970808087706;
	Matrix yPhi(42, 1);
	yPhi(1,1)=7101800.90695315;
	yPhi(2,1)=1293997.58115302;
	yPhi(3,1)=10114.014948955;
	yPhi(4,1)=573.068082065557;
	yPhi(5,1)=-3085.15736953138;
	yPhi(6,1)=-6736.03068347156;
	yPhi(7,1)=1.0000293469741;
	yPhi(8,1)=8.22733917593032e-06;
	yPhi(9,1)=2.17104932968693e-07;
	yPhi(10,1)=1.08925458231315e-05;
	yPhi(11,1)=3.04673932160225e-06;
	yPhi(12,1)=6.63504292706821e-08;
	yPhi(13,1)=8.22733944423959e-06;
	yPhi(14,1)=0.999986101965304;
	yPhi(15,1)=3.99927483270551e-08;
	yPhi(16,1)=3.04673960163327e-06;
	yPhi(17,1)=-5.1596062466179e-06;
	yPhi(18,1)=1.22075292404534e-08;
	yPhi(19,1)=2.17105640392839e-07;
	yPhi(20,1)=3.9992870847826e-08;
	yPhi(21,1)=0.999984551298692;
	yPhi(22,1)=6.63510875632706e-08;
	yPhi(23,1)=1.22076480274715e-08;
	yPhi(24,1)=-5.73276287738792e-06;
	yPhi(25,1)=5.38976081674752;
	yPhi(26,1)=1.47507305174403e-05;
	yPhi(27,1)=3.21241787851554e-07;
	yPhi(28,1)=1.00002936035846;
	yPhi(29,1)=8.19365458482084e-06;
	yPhi(30,1)=1.40504658112974e-07;
	yPhi(31,1)=1.47507306419397e-05;
	yPhi(32,1)=5.38968310056198;
	yPhi(33,1)=5.90697768748029e-08;
	yPhi(34,1)=8.19365482653896e-06;
	yPhi(35,1)=0.9999860891763;
	yPhi(36,1)=2.58022974647481e-08;
	yPhi(37,1)=3.21242427100724e-07;
	yPhi(38,1)=5.90698876854246e-08;
	yPhi(39,1)=5.38968032557769;
	yPhi(40,1)=1.4050537070756e-07;
	yPhi(41,1)=2.58024285760964e-08;
	yPhi(42,1)=0.999984550703337;

	Matrix expected(42, 1);
	expected(1,1)=573.068082065557;
	expected(2,1)=-3085.15736953138;
	expected(3,1)=-6736.03068347156;
	expected(4,1)=-7.53489822593659;
	expected(5,1)=-1.37294429126638;
	expected(6,1)=-0.0107597986473575;
	expected(7,1)=1.08925458231315e-05;
	expected(8,1)=3.04673932160225e-06;
	expected(9,1)=6.63504292706821e-08;
	expected(10,1)=2.02239897508587e-06;
	expected(11,1)=5.61811901849645e-07;
	expected(12,1)=4.39846387071934e-09;
	expected(13,1)=3.04673960163327e-06;
	expected(14,1)=-5.1596062466179e-06;
	expected(15,1)=1.22075292404534e-08;
	expected(16,1)=5.61812134084449e-07;
	expected(17,1)=-9.58613689243416e-07;
	expected(18,1)=8.05616500343474e-10;
	expected(19,1)=6.63510875632706e-08;
	expected(20,1)=1.22076480274715e-08;
	expected(21,1)=-5.73276287738792e-06;
	expected(22,1)=4.39895597958216e-09;
	expected(23,1)=8.0570607835305e-10;
	expected(24,1)=-1.06368693580442e-06;
	expected(25,1)=1.00002936035846;
	expected(26,1)=8.19365458482084e-06;
	expected(27,1)=1.40504658112974e-07;
	expected(28,1)=1.08999102436198e-05;
	expected(29,1)=3.02797128053784e-06;
	expected(30,1)=2.37068516291712e-08;
	expected(31,1)=8.19365482653896e-06;
	expected(32,1)=0.9999860891763;
	expected(33,1)=2.58022974647481e-08;
	expected(34,1)=3.02797160153579e-06;
	expected(35,1)=-5.16671243316801e-06;
	expected(36,1)=4.34211426867344e-09;
	expected(37,1)=1.4050537070756e-07;
	expected(38,1)=2.58024285760964e-08;
	expected(39,1)=0.999984550703337;
	expected(40,1)=2.37075280907946e-08;
	expected(41,1)=4.34223837651307e-09;
	expected(42,1)=-5.73302112206999e-06;

	Matrix result = VarEqn(x, yPhi);

	_assert(m_equals(result, expected, 1e-10));
	return 0;
}

int geodetic_01() {
	Matrix r(3);
	r(1)=34; r(2)=54; r(3)=3;
	auto [lon, lat, h] = Geodetic(r);

	double latexp = 1.569306928241210;
	double lonexp = 1.008874764538081;
	double hexp = -6.356748569071247e+06;
	_assert(fabs(lat-latexp) < 1e-10);
	_assert(fabs(h-hexp) < 1e-10);
	_assert(fabs(lon-lonexp) < 1e-10);
	return 0;
}

int angl_01() {
	Matrix v1(3);
	v1(1)=0; v1(2)=0; v1(3)=1;
	Matrix v2(3);
	v2(1)=1; v2(2)=0; v2(3)=0;
	double angle = angl(v1,v2);

	_assert(fabs(Constants::pi/2-angle) < 1e-10);
	return 0;
}

int elements_01() {
	Matrix V(6);
	V(1) = 6221397.62857869;
	V(2) = 2867713.77965738;
	V(3) = 3006155.98509949;

	V(4) = 4645.04725161806;
	V(5) = -2752.21591588204;
    V(6) = -7507.99940987031;

	auto [p, a, e, i, Omega, omega, M] = elements(V);
	double expectedP = 12001693.597214;
	double eA = 18943922.6607145;
	double eE = 0.605361104987026;
	double eI = 2.02656295535017;
	double eOmega = 3.35671076650829;
	double eomega = 2.73757289772562;
	double eM = 6.27144693341967;
	_assert(fabs(p-expectedP) < 1e-6);
	_assert(fabs(a-eA) < 1e-6);
	_assert(fabs(e-eE) < 1e-10);
	_assert(fabs(i-eI) < 1e-10);
	_assert(fabs(Omega-eOmega) < 1e-10);
	_assert(fabs(omega-eomega) < 1e-10);
	_assert(fabs(M-eM) < 1e-10);
	return 0;

}

int unit_01() {
	Matrix vec(3, 1);
	vec(1, 1)=20784239358.2383;
	vec(2, 1)=70646304145.8809;
	vec(3, 1)=-13038108211.875;
	Matrix exp(3);
	exp(1)=0.277917884211324; exp(2)=0.944651908456325; exp(3)=-0.174339959519679;

	Matrix result = unit(vec);
	_assert(m_equals(exp, result, 1e-10));

	return 0;
}

int gibbs_01() {
	Matrix r1(3); Matrix r2(3); Matrix r3(3);
	r1(1)=5720303.71012986;r1(2)=3152426.6965331;r1(3)=3750056.80416402;
	r2(1)=6221397.62857869;r2(2)=2867713.77965738;r2(3)=3006155.98509949;
	r3(1)=6699811.80976796;r3(2)=2569867.80763881;r3(3)=2154940.29542389;
	auto [v2, theta, theta1, copa, error] = gibbs(r1.transpose(), r2.transpose(), r3.transpose());
	Matrix ev2(3);
	ev2(1)=4645.04725161806; ev2(2)=-2752.21591588204; ev2(3)=-7507.99940987031;
	double etheta=0.125269502872995;
	double etheta1=0.136454013492468;
	double ecopa=0.00509723347775707;
	std::string eerror = "          ok";
	_assert(m_equals(v2, ev2.transpose(), 1e-9));
	_assert(fabs(theta-etheta) < 1e-10);
	_assert(fabs(theta1-etheta1) < 1e-10);
	_assert(fabs(copa-ecopa) < 1e-10);
	_assert(error==eerror);

	return 0;
}

int deinteg_01() {
	AuxParam.Mjd_UTC=4.974611128472211e+04;
	Matrix Y0_apr(6, 1);
	Y0_apr(1, 1) = 6.221397628578685e+06;
	Y0_apr(2, 1) = 2.867713779657379e+06;
	Y0_apr(3, 1) = 3.006155985099489e+06;
	Y0_apr(4, 1) = 4.645047251618060e+03;
	Y0_apr(5, 1) = -2.752215915882042e+03;
	Y0_apr(6, 1) = -7.507999409870306e+03;
	Matrix result = DEInteg(Accel, 0.0, -1.349999919533730e+02, 1e-13, 1e-6, 6, Y0_apr);

	Matrix expected(6, 1);
	expected(1, 1) = 5.542555937228607e+06;
	expected(2, 1) = 3.213514867349196e+06;
	expected(3, 1) = 3.990892975876853e+06;
	expected(4, 1) = 5.394068421663513e+03;
	expected(5, 1) = -2.365213378823415e+03;
	expected(6, 1) = -7.061845542002954e+03;
	// cout << result << endl;
	_assert(m_equals(result.transpose(), expected, 1e-9));
	return 0;

}

int hgibbs_01() {
	double Mjd1=4.974611015046295e+04;
	double Mjd2=4.974611128472211e+04;
	double Mjd3=4.974611253472231e+04;
	Matrix r1(3, 1);
	r1(1, 1)=5.720303710129856e+06;
	r1(2, 1)=3.152426696533103e+06;
	r1(3, 1)=3.750056804164019e+06;
	Matrix r2(3,1);
	r2(1, 1) = 6.221397628578685e+06;
	r2(2, 1) = 2.867713779657379e+06;
	r2(3, 1) = 3.006155985099489e+06;
	Matrix r3(3, 1);
	r3(1, 1) = 6.699811809767957e+06;
	r3(2, 1) = 2.569867807638814e+06;
	r3(3, 1) = 2.154940295423891e+06;

	Matrix ev2(3, 1);
	ev2(1, 1) = 4.796825169167805e+03;
	ev2(2, 1) = -2.839418128699730e+03;
	ev2(3, 1) = -7.741594338630217e+03;
	double etheta = 0.125269502872995;
	double etheta1 = 0.136454013492468;
	double ecopa = 0.005097233477757;
	std::string eerror = "   angl > 1ø";
	auto [v2, theta, theta1, copa, error] = hgibbs(r1, r2, r3, Mjd1, Mjd2, Mjd3);

	_assert(m_equals(v2, ev2, 1e-10));
	_assert(fabs(theta-etheta) < 1e-10);
	_assert(fabs(theta1-etheta1) < 1e-10);
	_assert(fabs(ecopa-copa) < 1e-10);
	_assert(error==eerror);

	return 0;
}

int anglesg_01() {
	double arg1=1.0559084894933;
	double arg2=1.36310214580757;
	double arg3=1.97615602688759;
	double arg4=0.282624656433946;
	double arg5=0.453434794338875;
	double arg6=0.586427138011591;
	double Mjd1=4.974611015046295e+04;
	double Mjd2=4.974611128472211e+04;
	double Mjd3=4.974611253472231e+04;
	Matrix Rs(3, 1);
	Rs(1, 1) = -5.512567840036068e+06;
	Rs(2, 1) = -2.196994446669333e+06;
	Rs(3, 1) = 2.330804966146887e+06;
	Matrix r2(3, 1);
    r2(1,1)=6221397.62857869;
    r2(2,1)=2867713.77965738;
    r2(3,1)=3006155.98509949;
    Matrix v2(3,1);
    v2(1,1)=4645.04725161806;
    v2(2,1)=-2752.21591588204;
    v2(3,1)=-7507.99940987031;

	auto [r2res, v2res] = anglesg(arg1, arg2, arg3, arg4, arg5, arg6, Mjd1, Mjd2, Mjd3, Rs, Rs, Rs);
	// cout << "r2res:" << r2res;
	// cout << "v2res:" << v2res;
	_assert(m_equals(r2, r2res, 1e-6));
	_assert(m_equals(v2, v2res, 1e-6));

	return 0;
}

int all_tests()
{
	eop19620101(21413); // c = 21413
	GGM03S();
	auxparam();
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
	_verify(cheb3d_02);
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
	_verify(ghamatrix_01);
	_verify(accel_01);
	_verify(vareqn_01);
	_verify(geodetic_01);
	_verify(angl_01);
	_verify(elements_01);
	_verify(unit_01);
	_verify(gibbs_01);
	_verify(deinteg_01);
	_verify(hgibbs_01);
	_verify(anglesg_01);

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
