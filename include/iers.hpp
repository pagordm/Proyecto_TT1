#ifndef _IERS_
#define _IERS_

#include "matrix.h"
#include <tuple>
#include <cmath>
#include "const.hpp"

std::tuple<double, double, double, double, double, double, double, double, double> IERS(Matrix& eop, double Mjd_UTC, char interp);

std::tuple<double, double, double, double, double, double, double, double, double> IERS(Matrix& eop, double Mjd_UTC);


#endif
