/**
 * @file gast.cpp
 * @author Pablo Gordillo Minchinela
 * @brief This file contains the implementation of the function gast.
 * @date 2025-05-20
 * 
 * 
 */
#include "..\include\gast.hpp"
/**
 * @brief Greenwich Apparent Sidereal Time
 * 
 * @param Mjd_UT1 Modified Julian Date UT1
 * @return double - GAST in [rad]
 */
double gast(double Mjd_UT1) {
    
    double gstime = fmod ( gmst(Mjd_UT1) + EqnEquinox(Mjd_UT1), 2*Constants::pi );
    if (gstime < 0) {
        gstime += 2*Constants::pi;
    }

    return gstime;
}
