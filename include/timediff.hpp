/**
 * @file timediff.hpp
 * @author Pablo Gordillo Minchinela
 * @brief This header contains the timediff function declaration.
 * @date 2025-05-20
 * 
 * 
 */
#ifndef _TIMEDIFF_
#define _TIMEDIFF_

#include <tuple>

std::tuple<double, double, double, double, double> timediff(double UT1_UTC, double TAI_UTC);

#endif
