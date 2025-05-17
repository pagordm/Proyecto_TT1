#include "..\include\measupdate.hpp"

std::tuple<Matrix&, Matrix&, Matrix&> MeasUpdate(Matrix x, double z, double g, double s, Matrix G, Matrix P, int n) {

    double m = 1;
    Matrix Inv_W = zeros(m,m);

    Inv_W(1,1) = s*s;    // Inverse weight (measurement covariance)
    

    // Kalman gain
    Matrix& K = P*G.transpose()*inv(Inv_W+G*P*G.transpose());

    // State update
    Matrix& newx = x + K*(z-g);

    // Covariance update
    Matrix& newP = (eye(n)-K*G)*P;

    return std::tie(K, newx, newP);

}
