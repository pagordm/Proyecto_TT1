#include "..\include\elements.hpp"

/**
 * @brief Computes the osculating Keplerian elements from the satellite state vector for elliptic orbits
 * 
 * This function calculates orbital elements from a given state vector. Note that this function
 * cannot be used with state vectors describing a circular or non-inclined orbit.
 * 
 * @param y State vector containing position and velocity (x,y,z,vx,vy,vz)
 * 
 * @return tuple containing:
 *         - double p: Semilatus rectum (meters)
 *         - double a: Semimajor axis
 *         - double e: Eccentricity
 *         - double i: Inclination (radians)
 *         - double Omega: Longitude of the ascending node (radians)
 *         - double omega: Argument of pericenter (radians)
 *         - double M: Mean anomaly (radians)
 */

std::tuple<double, double, double, double, double, double, double> elements(Matrix y) {
        
    double pi2, magh, H, p, Omega, i, u, R, a, eCosE, eSinE, e2, e, E, M, nu, omega;
    Matrix r, v, h;

    pi2 = Constants::pi2;
    
    r = y.extract_vector(1,3);                                        // Position
    v = y.extract_vector(4,6);                                        // Velocity
    
    h = cross(r,v);                                    // Areal velocity
    magh = norm(h);
    p = magh*magh/Constants::GM_Earth;
    H = norm(h);
    
    Omega = atan2 ( h(1), -h(2) );                     // Long. ascend. node 
    Omega = fmod(Omega,pi2);
    if (Omega < 0) {
        Omega += pi2;
    }
    i     = atan2 ( sqrt(h(1)*h(1)+h(2)*h(2)), h(3) ); // Inclination        
    u     = atan2 ( r(3)*H, -r(1)*h(2)+r(2)*h(1) );    // Arg. of latitude   
    
    R  = norm(r);                                      // Distance           
    
    a = 1.0/(2.0/R-dot(v,v)/Constants::GM_Earth);           // Semi-major axis    
    
    eCosE = 1.0-R/a;                                      // e*cos(E)           
    eSinE = dot(r,v)/sqrt(Constants::GM_Earth*a);       // e*sin(E)           
    
    e2 = eCosE*eCosE +eSinE*eSinE;
    e  = sqrt(e2);                                      // Eccentricity 
    E  = atan2(eSinE,eCosE);                            // Eccentric anomaly  
    
    M  = fmod(E-eSinE,pi2);                             // Mean anomaly
    if (M < 0) {
        M+=pi2;
    }
    
    nu = atan2(sqrt(1.0-e2)*eSinE, eCosE-e2);           // True anomaly
    
    omega = fmod(u-nu,pi2);                             // Arg. of perihelion 
    if (omega < 0) {
        omega += pi2;
    }
    
    return std::tie(p, a, e, i, Omega, omega, M);
    
}
