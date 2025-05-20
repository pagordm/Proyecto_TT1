/**
 * @file timeupdate.hpp
 * @author Pablo Gordillo Minchinela
 * @brief This header contains the timeupdate function declaration.
 * @date 2025-05-20
 * 
 * 
 */
#ifndef _TIMEUPDATE_
#define _TIMEUPDATE_

#include "matrix.h"

Matrix& TimeUpdate(Matrix& P, Matrix Phi, double Qdt);

Matrix& TimeUpdate(Matrix& P, Matrix Phi);

#endif
