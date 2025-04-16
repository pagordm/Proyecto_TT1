#include "..\include\sign_.hpp"

/**
 * @brief Returns absolute value of a with sign of b
 * 
 * @param a double to get absolute value of
 * @param b double to get sign of
 * @return double absolute value of a with sign of b
 */
double sign_(double a, double b) {
    a = a>0 ? a : -a; // absolute value of a
    if (b > 0) {
        return a;
    } else {
        return -a;
    }
}

