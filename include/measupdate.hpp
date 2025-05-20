/**
 * @file measupdate.hpp
 * @author Pablo Gordillo Minchinela
 * @brief This header contains the MeasUpdate function declaration.
 * @date 2025-05-20
 * 
 * 
 */
#ifndef _MEASUPDATE_
#define _MEASUPDATE_

#include "matrix.h"
#include <tuple>

std::tuple<Matrix&, Matrix&, Matrix&> MeasUpdate(Matrix x, double z, double g, double s, Matrix G, Matrix P, int n);

#endif