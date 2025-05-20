/**
 * @file R_z.cpp
 * @author Pablo Gordillo Minchinela
 * @brief This file contains the implementation of the function R_z.
 * @date 2025-05-20
 * 
 * 
 */
#include "..\include\R_z.hpp"
#include <cmath>
/**
 * @brief Generate a rotation matrix for a given angle around the z-axis
 * 
 * @param angle The angle of rotation in radians
 * @return Matrix& The rotation matrix
 */
Matrix& R_z(double angle) {
	double C = cos(angle);
	double S = sin(angle);
	Matrix& rotmat = zeros(3,3);

	rotmat(1,1) =      C;  rotmat(1,2) =   S;  rotmat(1,3) = 0.0;
	rotmat(2,1) = -1.0*S;  rotmat(2,2) =   C;  rotmat(2,3) = 0.0;
	rotmat(3,1) =    0.0;  rotmat(3,2) = 0.0;  rotmat(3,3) = 1.0;

	return (rotmat);
}
