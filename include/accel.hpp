/**
 * @file accel.hpp
 * @brief Declarations for satellite acceleration computations.
 *
 * This header file contains the declaration for the Accel function, which computes
 * the acceleration of an Earth-orbiting satellite, including the effects of Earth's
 * harmonic gravity field, gravitational perturbations from the Sun and Moon, solar
 * radiation pressure, and atmospheric drag.
 *
 * @author Pablo Gordillo
 */
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

Matrix& Accel(double x, Matrix Y);

#endif