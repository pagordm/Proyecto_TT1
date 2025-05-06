#include "..\include\matrix.h"
#include "..\include\global.hpp"
#include "..\include\legendre.hpp"
#include "..\include\JPL_Eph_DE430.hpp"
#include <iostream>

using namespace std;

int main() {
    eop19620101(4); // c = 21413
    GGM03S();
    DE430Coeff();

    auto [a,b,c,d,e,f,g,h,i,j,k] = JPL_Eph_DE430(49746.1);

    cout << "a:" << a << endl;
    cout << "b:" << b << endl;

    Matrix M1(3, 2);
	M1(1,1) = 5;
	
    Matrix M2(3, 2);
	M2(1,1) = -3;
	
    Matrix M3 = M1 - M2;

    Matrix V(3);
    V(2)=5;
    cout << "V\n" << V << "\n";
    cout << "M1\n" << M1 << "\n";
    cout << "M2\n" << M2 << "\n";
    cout << "M3\n" << M3 << "\n";
	
	cout << M1(1,1) << "\n";

    return 0;
}