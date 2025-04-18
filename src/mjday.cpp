#include "..\include\mjday.hpp"

/**
 * @brief Calculates Modified Julian Date from calendar date and time
 * 
 * @param year Year
 * @param month Month
 * @param day Day
 * @param hour Universal time hour
 * @param minute Universal time minute
 * @param second Universal time second
 * 
 * @return Modified Julian Date
 */

double Mjday(int year, int month, int day, int hour, int minute, double second) {


    double jd = 367.0 * year- floor( (7 * (year + floor( (month + 9) / 12.0) ) ) * 0.25 ) + floor( 275 * month / 9.0 ) + day + 1721013.5 + ( (second/60.0 + minute ) / 60.0 + hour ) / 24.0;

    double Mjd = jd-2400000.5;

    return Mjd;
}

double Mjday(int year, int month, int day) {
    return Mjday(year, month, day, 0, 0, 0);
}
