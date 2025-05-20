/**
 * @file nutmatrix.hpp
 * @author Pablo Gordillo Minchinela
 * @brief This header contains the NutMatrix function declaration.
 * @date 2025-05-20
 * 
 * 
 */
#ifndef _NUTMATRIX_
#define _NUTMATRIX_

#include "matrix.h"
#include "meanobliquity.hpp"
#include "nutangles.hpp"
#include "R_x.hpp"
#include "R_z.hpp"

Matrix& NutMatrix(double Mjd_TT);

#endif