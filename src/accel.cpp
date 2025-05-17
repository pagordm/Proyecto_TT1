#include "..\include\accel.hpp"

/**
 * @brief Computes the acceleration of an Earth orbiting satellite
 * 
 * Calculates acceleration due to:
 * - Earth's harmonic gravity field
 * - Gravitational perturbations of Sun and Moon 
 * - Solar radiation pressure
 * - Atmospheric drag
 *
 * @param x Terrestrial Time (Modified Julian Date)
 * @param Y Satellite state vector in the ICRF/EME2000 system
 * @return Matrix& Acceleration (a=d^2r/dt^2) in the ICRF/EME2000 system (COLUMN VECTOR)
 */
Matrix& Accel(double x, Matrix Y) {
    auto [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,AuxParam.Mjd_UTC + x/86400,'l');
    auto [UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
    long double Mjd_UT1 = AuxParam.Mjd_UTC + x/86400 + UT1_UTC/86400;
    long double Mjd_TT = AuxParam.Mjd_UTC + x/86400 + TT_UTC/86400;
    Matrix& P = PrecMatrix(Constants::MJD_J2000,Mjd_TT);
    Matrix& N = NutMatrix(Mjd_TT);
    Matrix& T = N * P;
    Matrix& E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;
    long double MJD_TDB = Mjday_TDB(Mjd_TT);
    auto [r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus,
    r_Neptune,r_Pluto,r_Moon,r_Sun] = JPL_Eph_DE430(MJD_TDB);
    bool transpose=false;
    if (Y.n_column==1) { //es un vector columna
        Y = Y.transpose(); //para evitar errores con extract_vector
        transpose=true;
    }
    // Acceleration due to harmonic gravity field
    Matrix &a = AccelHarmonic(Y.extract_vector(1, 3).transpose(), E, AuxParam.n, AuxParam.m);
    // Luni-solar perturbations
    if (AuxParam.sun) {
        a = a + AccelPointMass(Y.extract_vector(1, 3).transpose(),r_Sun,Constants::GM_Sun);
    }
    if (AuxParam.moon) {
        a = a + AccelPointMass(Y.extract_vector(1, 3).transpose(),r_Moon,Constants::GM_Moon);
    }
    // Planetary perturbations
    if (AuxParam.planets) {
        a = a + AccelPointMass(Y.extract_vector(1, 3).transpose(),r_Mercury,Constants::GM_Mercury);
        a = a + AccelPointMass(Y.extract_vector(1, 3).transpose(),r_Venus,Constants::GM_Venus);
        a = a + AccelPointMass(Y.extract_vector(1, 3).transpose(),r_Mars,Constants::GM_Mars);
        a = a + AccelPointMass(Y.extract_vector(1, 3).transpose(),r_Jupiter,Constants::GM_Jupiter);
        a = a + AccelPointMass(Y.extract_vector(1, 3).transpose(),r_Saturn,Constants::GM_Saturn);
        a = a + AccelPointMass(Y.extract_vector(1, 3).transpose(),r_Uranus,Constants::GM_Uranus);    
        a = a + AccelPointMass(Y.extract_vector(1, 3).transpose(),r_Neptune,Constants::GM_Neptune);
        a = a + AccelPointMass(Y.extract_vector(1, 3).transpose(),r_Pluto,Constants::GM_Pluto);
    }

    //dY = [Y(4:6);a];

    Matrix& dY = union_vector(Y.extract_vector(4, 6), a.transpose()).transpose();
    if (transpose) { 
        Y = Y.transpose(); //deshacemos el transpose
    }
    return dY;

}
