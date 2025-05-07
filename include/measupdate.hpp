#ifndef _MEASUPDATE_
#define _MEASUPDATE_

#include "matrix.h"
#include <tuple>

std::tuple<Matrix&, Matrix&, Matrix&> MeasUpdate(Matrix x, double z, double g, double s, Matrix G, Matrix P, int n);

#endif