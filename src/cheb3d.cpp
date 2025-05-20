/**
 * @file cheb3d.cpp
 * @author Pablo Gordillo Minchinela
 * @brief This file contains the implementation of the function Cheb3D
 * @date 2025-05-20
 * 
 * 
 */
#include "..\include\cheb3d.hpp"


/**
 * @brief Chebyshev approximation of 3-dimensional vectors
 *
 * @param t  Time value to evaluate at
 * @param N  Number of coefficients
 * @param Ta Begin interval
 * @param Tb End interval 
 * @param Cx Coefficients of Chebyshev polynomial (x-coordinate)
 * @param Cy Coefficients of Chebyshev polynomial (y-coordinate)
 * @param Cz Coefficients of Chebyshev polynomial (z-coordinate)
 * @return Matrix& Position vector at time t
 */
Matrix& cheb3d(double t, int N, double Ta, double Tb, Matrix& Cx, Matrix& Cy, Matrix& Cz) {
        
    if ( (t<Ta) || (Tb<t) ) {
        cout << "ERROR: Time out of range in Cheb3D::Value\n";
        exit(EXIT_FAILURE);
    }

    // Clenshaw algorithm
    double tau = (2*t-Ta-Tb)/(Tb-Ta);  
    Matrix f1 = zeros(1,3);
    Matrix f2 = zeros(1,3);
    Matrix old_f1 = zeros(1,3);

    for (int i = N; i >= 2; i--) {
        old_f1 = f1;
        Matrix coeffs(3);
        coeffs(1) = Cx(i);
        coeffs(2) = Cy(i);
        coeffs(3) = Cz(i);
        f1 = (f1*(2*tau))-f2+coeffs;
        f2 = old_f1;
    }
    Matrix aux(3);
    aux(1) = Cx(1);
    aux(2) = Cy(1);
    aux(3) = Cz(1);
    Matrix& ChebApp = (f1*tau)-f2+aux;
    return ChebApp;
}
