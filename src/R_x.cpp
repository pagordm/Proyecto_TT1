#include "..\include\R_x.hpp"
#include <cmath>
/**
 * @brief Generate a rotation matrix for a given angle around the x-axis
 * 
 * @param angle The angle of rotation in radians
 * @return Matrix& The rotation matrix
 */
Matrix& R_x(double angle) {
	double C = cos(angle);
	double S = sin(angle);
	Matrix& rotmat = zeros(3,3);

	rotmat(1,1) = 1.0;  rotmat(1,2) =    0.0;  rotmat(1,3) = 0.0;
	rotmat(2,1) = 0.0;  rotmat(2,2) =      C;  rotmat(2,3) =   S;
	rotmat(3,1) = 0.0;  rotmat(3,2) = -1.0*S;  rotmat(3,3) =   C;

	return (rotmat);
}
