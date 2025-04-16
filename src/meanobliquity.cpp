#include "..\include\meanobliquity.hpp"

/**
 * @brief Computes the mean obliquity of the ecliptic
 *
 * @param Mjd_TT Modified Julian Date (Terrestrial Time)
 * @return double - Mean obliquity of the ecliptic [rad]
 */
double MeanObliquity(double Mjd_TT) {
    double T = (Mjd_TT-Constants::MJD_J2000)/36525;

    double MOblq = Constants::Rad *( 84381.448/3600-(46.8150+(0.00059-0.001813*T)*T)*T/3600 );

    return MOblq;
}
