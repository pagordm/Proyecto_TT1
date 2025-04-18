#include "..\include\iers.hpp"
/**
 * @brief IERS: Management of IERS time and polar motion data
 * 
 * @param eop IERS Earth Orientation Parameters matrix
 * @param Mjd_UTC Modified Julian Date in UTC
 * @param interp Interpolation method ('l' for linear, 'n' for nearest)
 * @return std::tuple<double, double, double, double, double, double, double, double, double> 
 */
std::tuple<double, double, double, double, double, double, double, double, double> iers(Matrix& eop, double Mjd_UTC, char interp) {
    

    //return values
    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;


    if (interp =='l') {
        // linear interpolation
        double mjd = (floor(Mjd_UTC));


        //i = find(mjd==eop(4,:),1,'first');
        int i = 0;
        Matrix aux = eop.extract_row(4);
        for (int j = 1; j <= aux.n_column; j++) {
            if (aux(j) == mjd) {
                i = j;
                break;
            }
        }

        //preeop = eop(:,i);
        Matrix preeop = eop.extract_column(i);

        //nexteop = eop(:,i+1);
        Matrix nexteop = eop.extract_column(i+1);
        
        double mfme = 1440*(Mjd_UTC-floor(Mjd_UTC));
        double fixf = mfme/1440;
        // Setting of IERS Earth rotation parameters
        // (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
        x_pole  = preeop(5)+(nexteop(5)-preeop(5))*fixf;
        y_pole  = preeop(6)+(nexteop(6)-preeop(6))*fixf;
        UT1_UTC = preeop(7)+(nexteop(7)-preeop(7))*fixf;
        LOD     = preeop(8)+(nexteop(8)-preeop(8))*fixf;
        dpsi    = preeop(9)+(nexteop(9)-preeop(9))*fixf;
        deps    = preeop(10)+(nexteop(10)-preeop(10))*fixf;
        dx_pole = preeop(11)+(nexteop(11)-preeop(11))*fixf;
        dy_pole = preeop(12)+(nexteop(12)-preeop(12))*fixf;
        TAI_UTC = preeop(13);
        
        x_pole  = x_pole/Constants::Arcs;  // Pole coordinate [rad]
        y_pole  = y_pole/Constants::Arcs;  // Pole coordinate [rad]
        dpsi    = dpsi/Constants::Arcs;
        deps    = deps/Constants::Arcs;
        dx_pole = dx_pole/Constants::Arcs; // Pole coordinate [rad]
        dy_pole = dy_pole/Constants::Arcs; // Pole coordinate [rad]
    } else if (interp =='n') {
        double mjd = (floor(Mjd_UTC));
        //i = find(mjd==eop(4,:),1,'first');
        int i = 0;
        Matrix aux = eop.extract_row(4);
        for (int j = 1; j <= aux.n_column; j++) {
            if (aux(j) == mjd) {
                i = j;
                break;
            }
        }

        Matrix neweop = eop.extract_column(i);
        // Setting of IERS Earth rotation parameters
        // (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
        x_pole  = neweop(5)/Constants::Arcs;  // Pole coordinate [rad]
        y_pole  = neweop(6)/Constants::Arcs;  // Pole coordinate [rad]
        UT1_UTC = neweop(7);             // UT1-UTC time difference [s]
        LOD     = neweop(8);             // Length of day [s]
        dpsi    = neweop(9)/Constants::Arcs;
        deps    = neweop(10)/Constants::Arcs;
        dx_pole = neweop(11)/Constants::Arcs; // Pole coordinate [rad]
        dy_pole = neweop(12)/Constants::Arcs; // Pole coordinate [rad]
        TAI_UTC = neweop(13);            // TAI-UTC time difference [s]
    }

    return std::make_tuple(x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC);
}
/**
 * @brief IERS: Management of IERS time and polar motion data with nearest interpolation.
 * 
 * @param eop IERS Earth Orientation Parameters matrix
 * @param Mjd_UTC Modified Julian Date in UTC
 * @return std::tuple<double, double, double, double, double, double, double, double, double> 
 */
std::tuple<double, double, double, double, double, double, double, double, double> iers(Matrix& eop, double Mjd_UTC) {
    return iers(eop, Mjd_UTC, 'n');
}
