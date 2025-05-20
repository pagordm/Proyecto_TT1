/**
 * @file angl.cpp
 * @author Pablo Gordillo Minchinela
 * @brief This file contains the implementation of the function Angl.
 * @date 2025-05-20
 * 
 * 
 */
#include "..\include\angl.hpp"
/**
 * @brief Angle between two vectors
 * 
 * @param vec1 vector 1
 * @param vec2 vector 2
 * @return double - angle between the two vectors  -pi to pi
 */
double angl(Matrix vec1, Matrix vec2) {

    double small, undefined, magv1, magv2, temp, theta;

    small     = 0.00000001;
    undefined = 999999.1;

    magv1 = norm(vec1);
    magv2 = norm(vec2);

    if (magv1*magv2 > small^2) {
        temp= dot(vec1,vec2) / (magv1*magv2);
        if (abs( temp ) > 1.0) {
            double sign = temp > 0 ? 1 : -1;
            temp= sign* 1.0;
        }
        theta= acos( temp );
    } else {
        theta= undefined;
    }

    return theta;
}
