#include "..\include\position.hpp"

/**
 * @brief Calculate geocentric position from geodetic coordinates
 * 
 * @param lon Longitude in radians
 * @param lat Latitude in radians  
 * @param h Height in meters
 * @return Matrix& Position vector in Earth-fixed frame
 */
Matrix& position(double lon, double lat, double h) {
    // Earth parameters
    double R_equ = Constants::R_Earth;
    double f = Constants::f_Earth;

    // Calculate derived quantities
    double e2 = f * (2.0 - f);    // Square of eccentricity
    double CosLat = cos(lat);     // Cosine of geodetic latitude
    double SinLat = sin(lat);     // Sine of geodetic latitude

    // Calculate radius of curvature in prime vertical
    double N = R_equ / sqrt(1.0 - e2 * SinLat * SinLat);

    // Calculate position vector components
    Matrix& r = *(new Matrix(3));
    r(1) = (N + h) * CosLat * cos(lon);
    r(2) = (N + h) * CosLat * sin(lon);
    r(3) = ((1.0 - e2) * N + h) * SinLat;

    return r;
}
