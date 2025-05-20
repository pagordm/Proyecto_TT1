/**
 * @file g_accelharmonic.cpp
 * @author Pablo Gordillo Minchinela
 * @brief This file contains the implementation of the function G_AccelHarmonic.
 * @date 2025-05-20
 * 
 * 
 */
#include "..\include\g_accelharmonic.hpp"

/**
 * @brief Computes the gradient of the Earth's harmonic gravity field 
 * 
 * @param r Satellite position vector in the true-of-date system
 * @param U Transformation matrix to body-fixed system
 * @param n_max Gravity model degree
 * @param m_max Gravity model order
 * @return Matrix& - Gradient (G=da/dr) in the true-of-date system
 */
Matrix& G_AccelHarmonic(Matrix r, Matrix U, int n_max, int m_max) {
    if (r.n_row==1) {
        r = r.transpose(); //convertimos el vector fila en vector columna
    }
    double d = 1.0;   // Position increment [m]

    Matrix& G = zeros(3,3);
    Matrix dr = zeros(3,1);
    Matrix da;
    // Gradient
    for (int i=1; i <= 3; i++) {
        // Set offset in i-th component of the position vector
        dr = zeros(3,1);
        dr(i) = d;

        // Acceleration difference
        da = AccelHarmonic ( r+dr/2,U, n_max, m_max ) -
            AccelHarmonic ( r-dr/2,U, n_max, m_max );
        // Derivative with respect to i-th axis    
        G.assign_column(i, (da/d).transpose());
    }

    return G;
}
