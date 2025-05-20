/**
 * @file frac.cpp
 * @author Pablo Gordillo Minchinela
 * @brief This file contains the implementation of the function Frac.
 * @date 2025-05-20
 * 
 * 
 */
#include "..\include\frac.hpp"

/**
 * @brief Returns the fractional part of a number
 * 
 * @param x The number to get the fractional part of
 * @return double The fractional part of the number
 */
double Frac(double x) {
    return x - std::floor(x);
}


