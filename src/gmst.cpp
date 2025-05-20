/**
 * @file gmst.cpp
 * @author Pablo Gordillo Minchinela
 * @brief This file contains the implementation of the function gmst.
 * @date 2025-05-20
 * 
 * 
 */
#include "..\include\gmst.hpp"

/**
 * @brief Greenwich Mean Sidereal Time
 * 
 * @param Mjd_UT1 Modified Julian Date UT1
 * @return double - GMST in [rad]
 */
double gmst(double Mjd_UT1) {
    double Secs, MJD_J2000, UT1, Mjd_0, T_0, T, gmst, gmstime;
    Secs = 86400.0;                       // Seconds per day
    MJD_J2000 = 51544.5;

    Mjd_0 = floor(Mjd_UT1);
    UT1   = Secs*(Mjd_UT1-Mjd_0);         // [s]
    T_0   = (Mjd_0  -MJD_J2000)/36525.0;
    T     = (Mjd_UT1-MJD_J2000)/36525.0;

    gmst  = 24110.54841 + 8640184.812866*T_0 + 1.002737909350795*UT1
            + (0.093104-6.2e-6*T)*T*T;    // [s]

    gmstime = 2*Constants::pi*Frac(gmst/Secs);       // [rad], 0..2pi

    return gmstime;

}
