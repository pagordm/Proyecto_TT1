#include "..\include\hgibbs.hpp"

/**
 * @brief Implements the Herrick-Gibbs approximation for orbit determination
 * 
 * This function finds the middle velocity vector for 3 given position vectors using
 * the Herrick-Gibbs approximation method.
 * 
 * @param r1 First position vector in IJK coordinates (meters)
 * @param r2 Second position vector in IJK coordinates (meters) 
 * @param r3 Third position vector in IJK coordinates (meters)
 * @param Mjd1 Julian date of first sighting (days from 4713 BC)
 * @param Mjd2 Julian date of second sighting (days from 4713 BC) 
 * @param Mjd3 Julian date of third sighting (days from 4713 BC)
 * 
 * @return tuple containing:
 *         - Matrix& v2: Velocity vector for r2 in IJK coordinates (m/s)
 *         - double theta: Angle between vectors (radians)
 *         - double theta1: Angle between vectors (radians)
 *         - double copa: Angle between vectors (radians)
 *         - string error: Flag indicating success ("ok", etc)
 */
std::tuple<Matrix&, double, double, double, std::string> hgibbs(Matrix r1, Matrix r2, Matrix r3, double Mjd1, double Mjd2, double Mjd3) {
    
    double small, theta, theta1, copa, tover2, tolangle, dt21, dt31, dt32, magr1, magr2, magr3, term1, term2, term3;
    Matrix v2(3), p, q, w, pn, r1n;

    std::string error =  "          ok";
    theta = 0.0;
    theta1= 0.0;
    magr1 = norm( r1 );
    magr2 = norm( r2 );
    magr3 = norm( r3 );

    for (int i= 1; i<= 3; i++) {
        v2(i)= 0.0;
    }

    tolangle= 0.01745329251994;
    dt21= (Mjd2-Mjd1)*86400.0;
    dt31= (Mjd3-Mjd1)*86400.0;
    dt32= (Mjd3-Mjd2)*86400.0;

    p = cross( r2.transpose(),r3.transpose() ).transpose();
    pn = unit( p );
    r1n = unit( r1 );
    copa=  asin( dot( pn,r1n ) );

    if ( fabs( dot(r1n,pn) ) > 0.017452406 ) {
        error= "not coplanar";
    }

    theta  = angl( r1,r2 );
    theta1 = angl( r2,r3 );

    if ( (theta > tolangle) | (theta1 > tolangle) ) {
        error= "   angl > 1ø";
    }

    term1= -dt32*( 1.0/(dt21*dt31) + Constants::GM_Earth/(12.0*magr1*magr1*magr1) );
    term2= (dt32-dt21)*( 1.0/(dt21*dt32) + Constants::GM_Earth/(12.0*magr2*magr2*magr2) );
    term3=  dt21*( 1.0/(dt32*dt31) + Constants::GM_Earth/(12.0*magr3*magr3*magr3) );

    v2 =  r1*term1 + r2*term2 + r3*term3;

    return std::tie(v2, theta,theta1,copa, error);

}
