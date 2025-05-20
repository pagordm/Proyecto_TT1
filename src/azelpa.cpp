/**
 * @file azelpa.cpp
 * @author Pablo Gordillo Minchinela
 * @brief This file contains the implementation of the function AzElPa 
 * @date 2025-05-20
 * 
 * 
 */
#include "..\include\azelpa.hpp"

/**
 * @brief Computes azimuth, elevation and partials from local tangent coordinates
 *
 * @param s Topocentric local tangent coordinates (East-North-Zenith frame)
 * @return std::tuple containing:
 *         - double A: Azimuth [rad]
 *         - double E: Elevation [rad] 
 *         - Matrix dAds: Partials of azimuth w.r.t. s
 *         - Matrix dEds: Partials of elevation w.r.t. s
 */
std::tuple<double, double, Matrix&, Matrix&> AzElPa(Matrix& s) {
    double pi2 = Constants::pi2;

    double rho = sqrt(s(1)*s(1)+s(2)*s(2));

    // Angles
    double Az = atan2(s(1),s(2));

    if (Az<0.0) 
        Az = Az+pi2;
    

    double El = atan ( s(3) / rho );

    // Partials
    Matrix &dAds = zeros(3);
    dAds(1) = s(2)/(rho*rho);
    dAds(2) = -s(1)/(rho*rho);
    dAds(3) = 0.0;
    Matrix &dEds = zeros(3);
    dEds(1) = -s(1)*s(3)/rho;
    dEds(2) = -s(2)*s(3)/rho;
    dEds(3) = rho;
    dEds = dEds / dot(s,s);

    return std::tie(Az, El, dAds, dEds);
    
}
