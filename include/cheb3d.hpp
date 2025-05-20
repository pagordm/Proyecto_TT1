/**
 * @file cheb3d.hpp
 * @author Pablo Gordillo Minchinela
 * @brief This header contains the cheb3d function
 * 
 */
#ifndef _CHEB3D_
#define _CHEB3D_

#include "matrix.h"

Matrix& cheb3d(double t, int N, double Ta, double Tb, Matrix& Cx, Matrix& Cy, Matrix& Cz);

#endif

