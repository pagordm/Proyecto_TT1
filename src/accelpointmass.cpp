/**
 * @file accelpointmass.cpp
 * @author Pablo Gordillo Minchinela
 * @brief This file contains the implementation of the function AccelPointMass
 * @date 2025-05-20
 * 
 * 
 */
#include "..\include\accelpointmass.hpp"

/**
 * @brief Computes the perturbational acceleration due to a point mass
 * 
 * @param r Satellite position vector (column vector)
 * @param s Point mass position vector  (column vector)
 * @param GM Gravitational coefficient of point mass
 * @return Matrix& Acceleration (a=d^2r/dt^2) (column vector)
 */
Matrix& AccelPointMass(Matrix& r, Matrix& s, double GM) {
    // Relative position vector of satellite w.r.t. point mass
    Matrix d = r - s;

    // Acceleration
    
    Matrix& a = d / (pow(norm(d), 3.0)) + s / (pow(norm(s), 3.0));
    a = a * (-GM);
    return a;
}
