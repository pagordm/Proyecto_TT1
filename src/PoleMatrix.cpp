/**
 * @file PoleMatrix.cpp
 * @author Pablo Gordillo Minchinela
 * @brief This file contains the implementation of the function PoleMatrix.
 * @date 2025-05-20
 * 
 * 
 */
#include "..\include\PoleMatrix.hpp"

/**
 * @brief Transformation from pseudo Earth-fixed to Earth-fixed coordinates 
 * for a given date
 * 
 * @param xp Pole coordinte(xp,yp)
 * @param yp Pole coordinte(xp,yp)
 * @return Matrix& Pole matrix
 */
Matrix& PoleMatrix(double xp, double yp) {
    Matrix & polemat = R_y(-xp) * R_x(-yp);

    return polemat;
}
