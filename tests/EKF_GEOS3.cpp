#include "..\include\matrix.h"
#include "..\include\global.hpp"
#include "..\include\legendre.hpp"
#include "..\include\JPL_Eph_DE430.hpp"
#include "..\include\const.hpp"
#include "..\include\position.hpp"
#include "..\include\anglesg.hpp"
#include "..\include\accel.hpp"
#include "..\include\DEInteg.hpp"
#include "..\include\iers.hpp"
#include "..\include\timediff.hpp"
#include "..\include\vareqn.hpp"
#include "..\include\azelpa.hpp"
#include "..\include\timeupdate.hpp"
#include "..\include\measupdate.hpp"
#include <iostream>

using namespace std;

int main() {
    eop19620101(21413); // c = 21413
    GGM03S();
    DE430Coeff();
    auxparam();
    int nobs = 46;
    GEOS3(nobs);

    double sigma_range, sigma_az, sigma_el, lat, lon, alt, Mjd1, Mjd2, Mjd3, Mjd0, Mjd_UTC, Mjd_TT, Mjd_UT1, theta, Dist, Azim, Elev;
    int n_eqn, t, t_old;
    Matrix Rs, Y0_apr, Y, P, LT, yPhi, Phi, Y_old, U, r, s, dAdY, dEdY, dDds, dDdY, K, Y0, Y_true, dAds, dEds;
    sigma_range = 92.5;          // [m]
    sigma_az = 0.0224*Constants::Rad; // [rad]
    sigma_el = 0.0139*Constants::Rad; // [rad]

    // Kaena Point station
    lat = Constants::Rad*21.5748;     // [rad]
    lon = Constants::Rad*(-158.2706); // [rad]
    alt = 300.20;                // [m]
    Rs = position(lon, lat, alt).transpose();

    Mjd1 = obs(1,1);
    Mjd2 = obs(9,1);
    Mjd3 = obs(18,1);

    //auto [r2,v2] = anglesg(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),
    //                Mjd1,Mjd2,Mjd3,Rs,Rs,Rs);
    Matrix r2(3, 1);
    r2(1,1)=6221397.62857869;
    r2(2,1)=2867713.77965738;
    r2(3,1)=3006155.98509949;
    Matrix v2(3,1);
    v2(1,1)=4645.04725161806;
    v2(2,1)=-2752.21591588204;
    v2(3,1)=-7507.99940987031;


    // [r2,v2] = anglesdr(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),...
    //                    Mjd1,Mjd2,Mjd3,Rs,Rs,Rs);

    Y0_apr = union_vector(r2.transpose(), v2.transpose()).transpose();

    Mjd0 = Mjday(1995,1,29,02,38,0);

    AuxParam.Mjd_UTC = Mjd_UTC;

    Mjd_UTC = obs(9,1);


    n_eqn  = 6;
    Y = DEInteg(Accel,0,-(obs(9,1)-Mjd0)*86400.0,1e-13,1e-6,6,Y0_apr).transpose();   
    P = zeros(6, 6);
    
    for (int i=1; i <= 3; i++) {
        P(i,i)=1e8;
    }
    for (int i=4; i <= 6; i++) {
        P(i,i)=1e3;
    }
    LT = LTC(lon,lat);

    yPhi = zeros(42,1);
    Phi  = zeros(6, 6);

    // Measurement loop
    t = 0;

    for (int i=1; i <= nobs; i++) {    
        // Previous step
        t_old = t;
        Y_old = Y;
        
        // Time increment and propagation
        Mjd_UTC = obs(i,1);                       // Modified Julian Date
        t       = (Mjd_UTC-Mjd0)*86400.0;         // Time since epoch [s]
        auto [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,Mjd_UTC,'l');
        auto [UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
        Mjd_TT = Mjd_UTC + TT_UTC/86400.0;
        Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;
        AuxParam.Mjd_UTC = Mjd_UTC;
        AuxParam.Mjd_TT = Mjd_TT;
            
        for (int ii=1; ii <= 6; ii++) {
            yPhi(ii) = Y_old(ii);
            for (int j=1; j <= 6; j++) {  
                if (ii==j) {
                    yPhi(6*j+ii) = 1; 
                } else {
                    yPhi(6*j+ii) = 0;
                }
            }
        }
        yPhi = DEInteg(VarEqn,0,t-t_old,1e-13,1e-6,42,yPhi).transpose();
        //return 0;
        //cout << "yPhi: " << yPhi << endl;
        // Extract state transition matrices
        for (int j=1; j <= 6; j++) {
            //Phi(:,j) = yPhi(6*j+1:6*j+6);
            Phi.assign_column(j, yPhi.transpose().extract_vector(6*j+1, 6*j+6));
        }
        Y = DEInteg (Accel,0,t-t_old,1e-13,1e-6,6,Y_old).transpose();
        
        // Topocentric coordinates
        theta = gmst(Mjd_UT1);                    // Earth rotation
        U = R_z(theta);
        r = Y.transpose().extract_vector(1, 3).transpose(); //Y(1:3);
        s = LT*(U*r-Rs);                          // Topocentric position [m]
        // Time update
        P = TimeUpdate(P, Phi);
        // Azimuth and partials
        std::tie(Azim, Elev, dAds, dEds) = AzElPa(s);     // Azimuth, Elevation
        //dAdY = [dAds*LT*U,zeros(1,3)];
        dAdY = union_vector(dAds*LT*U, zeros(1,3));
        // Measurement update
        std::tie(K, Y, P) = MeasUpdate ( Y, obs(i,2), Azim, sigma_az, dAdY, P, 6 );
        // Elevation and partials
        r = Y.transpose().extract_vector(1, 3).transpose(); //Y(1:3);
        s = LT*(U*r-Rs);                          // Topocentric position [m]
        std::tie(Azim, Elev, dAds, dEds) = AzElPa(s);     // Azimuth, Elevation
        //dEdY = [dEds*LT*U,zeros(1,3)];
        dEdY = union_vector(dEds*LT*U, zeros(1,3));
        // Measurement update
        std::tie(K, Y, P) = MeasUpdate ( Y, obs(i,3), Elev, sigma_el, dEdY, P, 6 );

        // Range and partials
        r = Y.transpose().extract_vector(1, 3).transpose(); //Y(1:3);
        s = LT*(U*r-Rs);                          // Topocentric position [m]
        Dist = norm(s); dDds = (s/Dist).transpose();         // Range
        //dDdY = [dDds*LT*U,zeros(1,3)];
        dDdY = union_vector(dDds*LT*U,zeros(1,3));
        // Measurement update
        std::tie(K, Y, P) = MeasUpdate ( Y, obs(i,4), Dist, sigma_range, dDdY, P, 6 );
    }
    auto [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,obs(46,1),'l');
    auto [UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.Mjd_TT = Mjd_TT;
    cout << "Last DEInteg, Y: \n" << Y << endl;
    Y0 = DEInteg (Accel,0,-(obs(46,1)-obs(1,1))*86400.0,1e-13,1e-6,6,Y).transpose();
    
    //Y_true = [5753.173e3, 2673.361e3, 3440.304e3, 4.324207e3, -1.924299e3, -5.728216e3]';
    Y_true = zeros(6);
    Y_true(1) = 5753.173e3;
    Y_true(2) = 2673.361e3;
    Y_true(3) = 3440.304e3;
    Y_true(4) = 4.324207e3;
    Y_true(5) = -1.924299e3;
    Y_true(6) = -5.728216e3;
    Y_true = Y_true.transpose();
    
    printf("\nError of Position Estimation\n");
    printf("dX%10.1f [m]\n", Y0(1)-Y_true(1));
    printf("dY%10.1f [m]\n", Y0(2)-Y_true(2)); 
    printf("dZ%10.1f [m]\n", Y0(3)-Y_true(3));
    printf("\nError of Velocity Estimation\n");
    printf("dVx%8.1f [m/s]\n", Y0(4)-Y_true(4));
    printf("dVy%8.1f [m/s]\n", Y0(5)-Y_true(5));
    printf("dVz%8.1f [m/s]\n", Y0(6)-Y_true(6));


    return 0;
}