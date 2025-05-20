/**
 * @file accelharmonic.hpp
 * @author Pablo Gordillo Minchinela
 * @brief Harmonic acceleration functions
 * 
 */
#ifndef _ACCELHARMONIC_
#define _ACCELHARMONIC_

#include "matrix.h"
#include <cmath>

Matrix& AccelHarmonic(Matrix& r, Matrix& E, int n_max, int m_max);

#endif
