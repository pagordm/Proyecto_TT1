#include "..\include\geodetic.hpp"

/**
 * @brief Converts position vector to geodetic coordinates
 * 
 * @param r Position vector [m]
 * @return std::tuple<double,double,double> Tuple containing:
 *         - Longitude [rad]
 *         - Latitude [rad] 
 *         - Altitude [m]
 *
 * @note Last modified: 2015/08/12 M. Mahooti
 */
std::tuple<double, double, double> Geodetic(Matrix r) {
    
    double R_equ = Constants::R_Earth;
    double f     = Constants::f_Earth;

    double epsRequ, e2, X, Y, Z, rho2, lon, lat, h, dZ, ZdZ, Nh, SinPhi, N, dZ_new;

    epsRequ = Constants::eps*R_equ;        // Convergence criterion
    e2      = f*(2.0-f);        // Square of eccentricity

    X = r(1);                   // Cartesian coordinates
    Y = r(2);
    Z = r(3);
    rho2 = X*X + Y*Y;           // Square of distance from z-axis

    // Check validity of input data
    if (norm(r)==0.0) {
        printf ( " invalid input in Geodetic constructor\n" );
        
        lon = 0.0;
        lat = 0.0;
        h   = -R_equ;
    }

    // Iteration 
    dZ = e2*Z;

    while(1) {
        ZdZ    =  Z + dZ;
        Nh     =  sqrt ( rho2 + ZdZ*ZdZ ); 
        SinPhi =  ZdZ / Nh;                    // Sine of geodetic latitude
        N      =  R_equ / sqrt(1.0-e2*SinPhi*SinPhi);
        dZ_new =  N*e2*SinPhi;
        if ( fabs(dZ-dZ_new) < epsRequ ) {
            break;
        }
        dZ = dZ_new;
    }

    // Longitude, latitude, altitude
    lon = atan2 ( Y, X );
    lat = atan2 ( ZdZ, sqrt(rho2) );
    h   = Nh - N;

    return std::tie(lon, lat, h);

}
