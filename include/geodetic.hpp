/**
 * @file geodetic.hpp
 * @author Pablo Gordillo Minchinela
 * @brief This header contains the geodetic function declaration.
 * @date 2025-05-20
 * 
 * 
 */
#ifndef _GEODETIC_
#define _GEODETIC_

#include <tuple>
#include "matrix.h"
#include "const.hpp"

std::tuple<double, double, double> Geodetic(Matrix r);

#endif