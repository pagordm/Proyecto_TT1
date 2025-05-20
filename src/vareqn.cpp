/**
 * @file vareqn.cpp
 * @author Pablo Gordillo Minchinela
 * @brief This file contains the implementation of the function VarEqn.
 * @date 2025-05-20
 * 
 * 
 */
#include "..\include\vareqn.hpp"

/**
 * @brief Computes the variational equations, i.e. the derivative of the state vector and the state transition matrix
 * 
 * @param x Time since epoch in [s]
 * @param yPhi (6+36)-dim vector comprising the state vector (y) and the state transition matrix (Phi) in column wise storage order
 * @return Matrix& Reference to yPhip, the derivative of yPhi
 * 
 * @note Last modified: 2015/08/12 M. Mahooti
 */
Matrix& VarEqn(double x, Matrix yPhi) {
    if (yPhi.n_row==1) { // Si es un vector fila lo trasponemos a un vector columna
        yPhi=yPhi.transpose();
    }
    auto [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,AuxParam.Mjd_UTC,'l');
    auto [UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
    double Mjd_UT1 = AuxParam.Mjd_TT + (UT1_UTC-TT_UTC)/86400;

    Matrix P, N, T, E, r, v, Phi, a, G, dfdy, Phip;

    // Transformation matrix
    P = PrecMatrix(Constants::MJD_J2000,AuxParam.Mjd_TT + x/86400);
    N = NutMatrix(AuxParam.Mjd_TT + x/86400);
    T = N * P;
    E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    // State vector components
    r = yPhi.extract_vector(1,3);
    v = yPhi.extract_vector(4,6);
    Phi = zeros(6, 6);

    // State transition matrix
    for (int j=1; j <= 6; j++) {
        //Phi(:,j) = yPhi(6*j+1:6*j+6);
        Phi.assign_column(j, yPhi.extract_vector(6*j+1, 6*j+6));
    }

    // Acceleration and gradient
    a = AccelHarmonic ( r, E, AuxParam.n, AuxParam.m ).transpose();
    G = G_AccelHarmonic ( r.transpose(), E, AuxParam.n, AuxParam.m ); //r.transpose() needed?

    // Time derivative of state transition matrix
    Matrix& yPhip = zeros(42,1);
    dfdy = zeros(6, 6);

    for (int i=1; i <= 3; i++) {
        for (int j=1; j <= 3; j++) {
            dfdy(i,j) = 0.0;                 // dv/dr(i,j)
            dfdy(i+3,j) = G(i,j);            // da/dr(i,j)
            if ( i==j ) {
                dfdy(i,j+3) = 1;
            } else {
                dfdy(i,j+3) = 0;             // dv/dv(i,j)
            }
            dfdy(i+3,j+3) = 0.0;             // da/dv(i,j)
        }
    }

    Phip = dfdy*Phi;

    // Derivative of combined state vector and state transition matrix
    for (int i=1; i <= 3; i++) {
        yPhip(i)   = v(i);                 // dr/dt(i)
        yPhip(i+3) = a(i);                 // dv/dt(i)
    }

    for (int i=1; i <= 6; i++) {
        for (int j=1; j<=6; j++) {
            yPhip(6*j+i) = Phip(i,j);     // dPhi/dt(i,j)
        }
    }
    
    return yPhip;

}
