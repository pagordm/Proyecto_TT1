/**
 * @file measupdate.cpp
 * @author Pablo Gordillo Minchinela
 * @brief This file contains the implementation of the function MeasUpdate
 * @date 2025-05-20
 * 
 * 
 */
#include "..\include\measupdate.hpp"
/**
 * @brief Measurement update
 * 
 * @param x a matrix
 * @param z double
 * @param g double
 * @param s double
 * @param G a matrix
 * @param P a matrix
 * @param n integer
 * @return std::tuple<Matrix&, Matrix&, Matrix&> - K, new X, new P
 */
std::tuple<Matrix&, Matrix&, Matrix&> MeasUpdate(Matrix x, double z, double g, double s, Matrix G, Matrix P, int n) {
    // cout << "MeasUpdate inputs:\n x:\n" << x << "z: " << z << "\ng: " << g << "\ns: " << s << "\nG:\n" << G << "P:\n" << P << "n: " << n << endl;
    if(x.n_row==1){
		x=x.transpose();
	}
    double m = 1;
    Matrix Inv_W(1);

    Inv_W(1) = s*s;    // Inverse weight (measurement covariance)
    

    // Kalman gain
    Matrix& K = P*G.transpose()*inv(Inv_W+G*P*G.transpose());

    // State update
    Matrix& newx = x + K*(z-g);

    // Covariance update
    Matrix& newP = (eye(n)-K*G)*P;
    // cout << "MeasUpdate outputs:\n K:\n" << K << "newx:\n" <<newx << "newP:\n" << newP << endl;
    return std::tie(K, newx, newP);

}
