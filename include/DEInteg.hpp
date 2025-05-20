/**
 * @file DEInteg.hpp
 * @author Pablo Gordillo Minchinela
 * @brief This function contains the DEInteg integration function.
 * @date 2025-05-20
 * 
 * 
 */
#ifndef _DEINTEG_
#define _DEINTEG_

#include "matrix.h"
#include "const.hpp"
#include "sign_.hpp"
#include <cmath>

Matrix& DEInteg(Matrix& f(double t, Matrix y), double t, double tout, double relerr, double abserr, int n_eqn, Matrix &y);

#endif