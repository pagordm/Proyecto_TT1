#include "..\include\LTC.hpp"

/**
 * @brief Transformation from Greenwich meridian system to 
 * local tangent coordinates
 * 
 * @param lon Geodetic East longitude [rad]
 * @param lat Geodetic latitude [rad]
 * @return Matrix& Rotation matrix from the Earth equator and Greenwich meridian
 * to the local tangent (East-North-Zenith) coordinate system
 */
Matrix& LTC(double lon, double lat) {
    Matrix& M = R_y(-1.0*lat)*R_z(lon);
    double Aux;
    for (int j=1; j <= 3; j++) {
        Aux=M(1,j); M(1,j)=M(2,j); M(2,j)=M(3,j); M(3,j)= Aux;
    }

    return M;
}
