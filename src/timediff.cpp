/**
 * @file timediff.cpp
 * @author Pablo Gordillo Minchinela
 * @brief This file contains the implementation of the function timediff.
 * @date 2025-05-20
 * 
 * 
 */
#include "..\include\timediff.hpp"
/**
 * @brief Calculates time differences between different time systems
 * 
 * @param UT1_UTC UT1-UTC time difference [s]
 * @param TAI_UTC TAI-UTC time difference [s]
 * @return std::tuple<double, double, double, double, double> 
 *         - UT1_TAI: UT1-TAI time difference [s]
 *         - UTC_GPS: UTC-GPS time difference [s]
 *         - UT1_GPS: UT1-GPS time difference [s]
 *         - TT_UTC: TT-UTC time difference [s]
 *         - GPS_UTC: GPS-UTC time difference [s]
 */
std::tuple<double, double, double, double, double> timediff(double UT1_UTC, double TAI_UTC) {
    double TT_TAI  = +32.184;          // TT-TAI time difference [s]
    double GPS_TAI = -19.0;            // GPS-TAI time difference [s]
    //Estos dos no se usan
    //double TT_GPS  =  TT_TAI-GPS_TAI;  // TT-GPS time difference [s]
    //double TAI_GPS = -GPS_TAI;         // TAI-GPS time difference [s]
    double UT1_TAI = UT1_UTC-TAI_UTC;  // UT1-TAI time difference [s]
    double UTC_TAI = -TAI_UTC;         // UTC-TAI time difference [s]
    double UTC_GPS = UTC_TAI-GPS_TAI;  // UTC_GPS time difference [s]
    double UT1_GPS = UT1_TAI-GPS_TAI;  // UT1-GPS time difference [s]
    double TT_UTC  = TT_TAI-UTC_TAI;   //  TT-UTC time difference [s]
    double GPS_UTC = GPS_TAI-UTC_TAI;  // GPS-UTC time difference [s]
    
    return std::tie(UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);

}


