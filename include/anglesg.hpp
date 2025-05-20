/**
 * @file anglesg.hpp
 * @author Pablo Gordillo Minchinela
 * @brief Anglesg function for EKF
 * 
 */
#ifndef _ANGLESG_
#define _ANGLESG_

#include <tuple>
#include "matrix.h"
#include "const.hpp"
#include "global.hpp"
#include "geodetic.hpp"
#include "LTC.hpp"
#include "iers.hpp"
#include "timediff.hpp"
#include "PrecMatrix.hpp"
#include "nutmatrix.hpp"
#include "PoleMatrix.hpp"
#include "GHAMatrix.hpp"
#include "angl.hpp"
#include "gibbs.hpp"
#include "hgibbs.hpp"
#include "elements.hpp"
#include "rpoly.h"

std::tuple<Matrix&, Matrix&> anglesg(double az1, double az2, double az3, double el1, double el2, double el3, double Mjd1, double Mjd2, double Mjd3, Matrix Rs1, Matrix Rs2, Matrix Rs3);

#endif