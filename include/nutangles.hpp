/**
 * @file nutangles.hpp
 * @author Pablo Gordillo Minchinela
 * @brief This header contains the NutAngles function declaration.
 * @date 2025-05-20
 * 
 * 
 */
#ifndef _NUTANGLES_
#define _NUTANGLES_

#include <tuple>
#include "const.hpp"
#include "matrix.h"
#include <cmath>


std::tuple<double, double> NutAngles(double Mjd_TT);

#endif