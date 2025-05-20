/**
 * @file position.hpp
 * @author Pablo Gordillo Minchinela
 * @brief This header contains the position function declaration.
 * @date 2025-05-20
 * 
 * 
 */
#ifndef _POSITION_
#define _POSITION_

#include "..\include\matrix.h"
#include "..\include\const.hpp"
#include <cmath>

using namespace std;

Matrix& position(double lon, double lat, double h);

#endif