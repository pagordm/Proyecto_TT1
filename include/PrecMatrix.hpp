/**
 * @file PrecMatrix.hpp
 * @author Pablo Gordillo Minchinela
 * @brief This header contains the PrecMatrix function declaration.
 * @date 2025-05-20
 * 
 * 
 */
#ifndef _PRECMATRIX_
#define _PRECMATRIX_

#include "matrix.h"
#include "const.hpp"
#include "R_z.hpp"
#include "R_y.hpp"

Matrix& PrecMatrix(double Mjd_1, double Mjd_2);

#endif