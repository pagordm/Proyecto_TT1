/**
 * @file legendre.hpp
 * @author Pablo Gordillo Minchinela
 * @brief This header contains the Legendre function declaration.
 * @date 2025-05-20
 * 
 * 
 */
#ifndef _LEGENDRE_
#define _LEGENDRE_

#include "matrix.h"
#include <tuple>

std::tuple<Matrix&, Matrix&> Legendre(int n, int m, double fi);

#endif

