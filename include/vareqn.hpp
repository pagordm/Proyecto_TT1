#ifndef _VAREQN_
#define _VAREQN_

#include "matrix.h"
#include "iers.hpp"
#include "timediff.hpp"
#include "PrecMatrix.hpp"
#include "nutmatrix.hpp"
#include "PoleMatrix.hpp"
#include "GHAMatrix.hpp"
#include "accelharmonic.hpp"
#include "g_accelharmonic.hpp"
#include "global.hpp"
#include "const.hpp"

Matrix& VarEqn(double x, Matrix& yPhi);

#endif