/**
 * @file sign_.cpp
 * @author Pablo Gordillo Minchinela
 * @brief This file contains the implementation of the function sign_.
 * @date 2025-05-20
 * 
 * 
 */
#include "..\include\sign_.hpp"

/**
 * @brief Returns absolute value of a with sign of b
 * 
 * @param a double to get absolute value of
 * @param b double to get sign of
 * @return double absolute value of a with sign of b
 */
double sign_(double a, double b) {
    if (b > 0) {
        return fabs(a);
    } else {
        return -fabs(a);
    }
}

