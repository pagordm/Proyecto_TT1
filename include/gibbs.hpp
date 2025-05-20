/**
 * @file gibbs.hpp
 * @author Pablo Gordillo Minchinela
 * @brief This header contains the gibbs function declaration.
 * @date 2025-05-20
 * 
 * 
 */
#ifndef _GIBBS_
#define _GIBBS_

#include <tuple>
#include "matrix.h"
#include "angl.hpp"
#include "unit.hpp"
#include "const.hpp"
#include <cmath>

std::tuple<Matrix&, double, double, double, std::string> gibbs(Matrix r1, Matrix r2, Matrix r3);

#endif