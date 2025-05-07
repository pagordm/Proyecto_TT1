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


