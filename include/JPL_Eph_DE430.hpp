/**
 * @file JPL_Eph_DE430.hpp
 * @author Pablo Gordillo Minchinela
 * @brief This header contains the JPL_Eph_DE430 function declaration.
 * @date 2025-05-20
 * 
 * 
 */
#ifndef _JPL_
#define _JPL_

#include <tuple>
#include "matrix.h"
#include "global.hpp"
#include "cheb3d.hpp"

std::tuple<Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&> JPL_Eph_DE430(double Mjd_TDB);


#endif