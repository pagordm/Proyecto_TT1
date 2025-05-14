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