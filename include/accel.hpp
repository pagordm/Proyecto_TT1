#ifndef _ACCEL_
#define _ACCEL_

#include "matrix.h"
#include "iers.hpp"
#include "timediff.hpp"
#include "PrecMatrix.hpp"
#include "nutmatrix.hpp"
#include "PoleMatrix.hpp"
#include "GHAMatrix.hpp"
#include "mjday_tdb.hpp"
#include "JPL_Eph_DE430.hpp"
#include "accelharmonic.hpp"
#include "accelpointmass.hpp"
#include "global.hpp"
#include "const.hpp"

Matrix& Accel(double x, Matrix& Y);

#endif