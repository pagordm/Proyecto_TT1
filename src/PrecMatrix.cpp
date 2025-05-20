/**
 * @file PrecMatrix.cpp
 * @author Pablo Gordillo Minchinela
 * @brief This file contains the implementation of the function PrecMatrix.
 * @date 2025-05-20
 * 
 * 
 */
#include "..\include\PrecMatrix.hpp"

Matrix& PrecMatrix(double Mjd_1, double Mjd_2) {
    double T  = (Mjd_1-Constants::MJD_J2000)/36525;
    double dT = (Mjd_2-Mjd_1)/36525;

    // Precession angles
    double zeta  =  ( (2306.2181+(1.39656-0.000139*T)*T)+
                ((0.30188-0.000344*T)+0.017998*dT)*dT )*dT/Constants::Arcs;
    double z     =  zeta + ( (0.79280+0.000411*T)+0.000205*dT)*dT*dT/Constants::Arcs;
    double theta =  ( (2004.3109-(0.85330+0.000217*T)*T)-
                ((0.42665+0.000217*T)+0.041833*dT)*dT )*dT/Constants::Arcs;

    // Precession matrix
    Matrix& PrecMat = R_z(-z) * R_y(theta) * R_z(-zeta);

    return PrecMat;
}
