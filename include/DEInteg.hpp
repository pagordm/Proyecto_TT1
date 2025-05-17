#ifndef _DEINTEG_
#define _DEINTEG_

#include "matrix.h"
#include "const.hpp"
#include "sign_.hpp"
#include <cmath>

Matrix& DEInteg(Matrix& f(double t, Matrix y), double t, double tout, double relerr, double abserr, int n_eqn, Matrix &y);

#endif