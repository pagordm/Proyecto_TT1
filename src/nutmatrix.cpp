/**
 * @file nutmatrix.cpp
 * @author Pablo Gordillo Minchinela
 * @brief This file contains the implementation of the function NutMatrix.
 * @date 2025-05-20
 * 
 * 
 */
#include "..\include\nutmatrix.hpp"

/**
 * @brief Transformation from mean to true equator and equinox
 * 
 * @param Mjd_TT Modified Julian Date (Terrestrial Time)
 * @return Matrix& Nutation matrix
 */
Matrix& NutMatrix(double Mjd_TT) {
    
    // Mean obliquity of the ecliptic
    double eps = MeanObliquity (Mjd_TT);

    // Nutation in longitude and obliquity
    auto [dpsi, deps] = NutAngles (Mjd_TT);

    // Transformation from mean to true equator and equinox
    Matrix& NutMat = R_x(-eps-deps)*R_z(-dpsi)*R_x(+eps);

    return NutMat;
}
