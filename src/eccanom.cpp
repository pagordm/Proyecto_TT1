#include "../include/eccanom.hpp"

#include <cmath>


/**
 * @brief Computes the eccentric anomaly for elliptic orbits
 *
 * @param M Mean anomaly in [rad]
 * @param e Eccentricity of the orbit [0,1]
 * @return Eccentric anomaly in [rad]
 */
double EccAnom(double M, double e) {
    int maxit = 15;
    int i = 1;
    // Starting value
    M = std::fmod(M, 2.0*Constants::pi);
    if (M < 0) {
        M += 2.0*Constants::pi;
    }
    double E = 0;
    if (e<0.8)
        E = M; 
    else
        E = Constants::pi;

    double f = E - e*sin(E) - M;
    E = E - f / ( 1.0 - e*cos(E) );

    // Iteration
    while (fabs(f) > 1e2*1e-10) {//1e2*Constants::eps) {  
        f = E - e*sin(E) - M;
        E = E - f / ( 1.0 - e*cos(E) );
        i = i+1;
        if (i==maxit) {
            std::cout << "convergence problems in EccAnom" << std::endl;
            exit(EXIT_FAILURE);
        }  
    }

    return E;
}
