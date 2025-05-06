#include "..\include\eqnequinox.hpp"

/**
 * @brief Computation of the equation of the equinoxes
 * 
 * @param Mjd_TT Modified Julian Date (Terrestrial Time)
 * @return double Equation of the equinoxes
 * 
 * @note The equation of the equinoxes dpsi*cos(eps) is the right ascension of
 *  the mean equinox referred to the true equator and equinox and is equal
 *  to the difference between apparent and mean sidereal time.
 */
double EqnEquinox(double Mjd_TT) {

    
    // Nutation in longitude and obliquity
    auto [dpsi, deps] = NutAngles (Mjd_TT);

    // Equation of the equinoxes
    double EqE = dpsi * cos ( MeanObliquity(Mjd_TT) );

    return EqE;

}