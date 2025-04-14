#include "..\include\accelpointmass.hpp"

/**
 * @brief Computes the perturbational acceleration due to a point mass
 * 
 * @param r Satellite position vector
 * @param s Point mass position vector  
 * @param GM Gravitational coefficient of point mass
 * @return Matrix& Acceleration (a=d^2r/dt^2)
 */
Matrix& accelpointmass(Matrix& r, Matrix& s, double GM) {
    // Relative position vector of satellite w.r.t. point mass
    Matrix& d = r - s;

    // Acceleration
    Matrix& a = *(new Matrix(3));
    a = (d / pow(norm(d), 3) + s / pow(norm(s), 3)) * (-GM);

    return a;
}
