/**
 * @file hgibbs.hpp
 * @author Pablo Gordillo Minchinela
 * @brief This header contains the hgibbs function declaration.
 * @date 2025-05-20
 * 
 * 
 */
#ifndef _HGIBBS_
#define _HGIBBS_

#include <tuple>
#include "matrix.h"
#include "angl.hpp"
#include "unit.hpp"
#include "const.hpp"
#include <cmath>

std::tuple<Matrix&, double, double, double, std::string> hgibbs(Matrix r1, Matrix r2, Matrix r3, double Mjd1, double Mjd2, double Mjd3);

#endif