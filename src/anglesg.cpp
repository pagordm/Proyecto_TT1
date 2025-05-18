#include "..\include\anglesg.hpp"

/**
 * @brief Solves the problem of orbit determination using three optical sightings
 * 
 * @param az1 Azimuth at t1 (rad)
 * @param az2 Azimuth at t2 (rad) 
 * @param az3 Azimuth at t3 (rad)
 * @param el1 Elevation at t1 (rad)
 * @param el2 Elevation at t2 (rad)
 * @param el3 Elevation at t3 (rad)
 * @param Mjd1 Modified julian date of t1
 * @param Mjd2 Modified julian date of t2 
 * @param Mjd3 Modified julian date of t3
 * @param Rs1 IJK site1 position vector (m)
 * @param Rs2 IJK site2 position vector (m)
 * @param Rs3 IJK site3 position vector (m)
 * @return std::tuple<Matrix&,Matrix&> Position vector at t2 (m), Velocity vector at t2 (m/s)
 * 
 * @note Last modified: 2015/08/12 M. Mahooti
 */

std::tuple<Matrix&, Matrix&> anglesg(double az1, double az2, double az3, double el1, double el2, double el3, double Mjd1, double Mjd2, double Mjd3, Matrix Rs1, Matrix Rs2, Matrix Rs3) {

    Matrix L1(3); Matrix L2(3); Matrix L3(3);
    Matrix M1, M2, M3, Lb1, Lb2, Lb3;
    Matrix Lm1(3), Lm2(3), Lm3(3);
    double Mjd_TT, Mjd_UT1, Mjd_UTC;
    Matrix P, N, T, E;
    double tau1, tau3;
    double a1, a3;
    double b1, b3;
    Matrix D;
    double d1s, d2s;
    double Ccye;
    Matrix rootarr;
    double bigr2;
    double u;
    double C1, C2, C3;
    Matrix temp;
    double rho1, rho2, rho3;
    double rhoold1, rhoold2, rhoold3;
    int ll;
    Matrix r1(3), r2(3), r3(3);
    double magr1, magr2, magr3;
    Matrix v2(3);
    double theta, theta1, copa;
    std::string error;
    double p, a, e, i, Omega, omega, M;
    double f1, f3, g1, g3;
    double rdot, udot;
    double H1, H2, H3;
    double G1, G2, G3;
    double D1, D2, D3;
    double tausqr;
    double UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC;
    double x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC;
    // cout << "1" << endl;
    L1(1)=cos(el1)*sin(az1); L1(2)= cos(el1)*cos(az1); L1(3)= sin(el1);;
    L2(1)=cos(el2)*sin(az2); L2(2)= cos(el2)*cos(az2); L2(3)= sin(el2);;
    L3(1)=cos(el3)*sin(az3); L3(2)= cos(el3)*cos(az3); L3(3)= sin(el3);;
    // cout << "2" << endl;
    auto [lon1, lat1, h1] = Geodetic(Rs1);
    auto [lon2, lat2, h2] = Geodetic(Rs2);
    auto [lon3, lat3, h3] = Geodetic(Rs3);
    // cout << "lon1: "<< lon1 << ", lat1: " << lat1 << endl;
    M1 = LTC(lon1, lat1);
    M2 = LTC(lon2, lat2);
    M3 = LTC(lon3, lat3);
    // cout << "m1: \n" << M1;
    

    // cout << "4" << endl;
    // body-fixed system
    Lb1 = M1.transpose()*L1.transpose();
    Lb2 = M1.transpose()*L2.transpose();
    Lb3 = M1.transpose()*L3.transpose();
    // cout << "5" << endl;
    // mean of date system (J2000)
    
    Mjd_UTC = Mjd1;
    std::tie(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC) = IERS(eopdata,Mjd_UTC,'l');
    std::tie(UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC) = timediff(UT1_UTC,TAI_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;
    // cout << "6" << endl;
    P = PrecMatrix(Constants::MJD_J2000,Mjd_TT);
    N = NutMatrix(Mjd_TT);
    T = N * P;
    E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;
    // cout << "7" << endl;
    Lm1 = E.transpose()*Lb1;
    Rs1 = E.transpose()*Rs1;
    // cout << "8" << endl;
    
    Mjd_UTC = Mjd2;
    std::tie(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC) = IERS(eopdata,Mjd_UTC,'l');
    std::tie(UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC) = timediff(UT1_UTC,TAI_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400;

    P = PrecMatrix(Constants::MJD_J2000,Mjd_TT);
    N = NutMatrix(Mjd_TT);
    T = N * P;
    E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    Lm2 = E.transpose()*Lb2;
    Rs2 = E.transpose()*Rs2;
    
    Mjd_UTC = Mjd3;
    std::tie(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC) = IERS(eopdata,Mjd_UTC,'l');
    std::tie(UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC) = timediff(UT1_UTC,TAI_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;

    P = PrecMatrix(Constants::MJD_J2000,Mjd_TT);
    N = NutMatrix(Mjd_TT);
    T = N * P;
    E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    Lm3 = E.transpose()*Lb3;
    Rs3 = E.transpose()*Rs3;

    // geocentric inertial position
    tau1 = (Mjd1-Mjd2)*86400;
    tau3 = (Mjd3-Mjd2)*86400;

    a1 = tau3/(tau3-tau1);
    a3 =-tau1/(tau3-tau1);

    b1 = tau3/(6*(tau3-tau1))*(pow((tau3-tau1),2.0)-pow(tau3,2.0));
    b3 =-tau1/(6*(tau3-tau1))*(pow((tau3-tau1),2.0)-pow(tau1,2.0));

    Matrix aux1(3,3), aux2(3,3);
    aux1.assign_column(1, Lm1.transpose());
    aux1.assign_column(2, Lm2.transpose());
    aux1.assign_column(3, Lm3.transpose());
    
    aux2.assign_column(1, Rs1.transpose());
    aux2.assign_column(2, Rs2.transpose());
    aux2.assign_column(3, Rs2.transpose());

    D = inv(aux1)*aux2;

    d1s = D(2,1)*a1-D(2,2)+D(2,3)*a3;
    d2s = D(2,1)*b1+D(2,3)*b3;

    Ccye = 2*dot(Lm2,Rs2);
    
    double poly[10] = {
        1.0,  // R2^8... polynomial
        0.0,
        -(pow(d1s,2) + d1s*Ccye + pow(norm(Rs2),2)),
        0.0,
        0.0,
        -Constants::GM_Earth*(d2s*Ccye + 2*d1s*d2s),
        0.0,
        0.0,
        -pow(Constants::GM_Earth,2)*pow(d2s,2)
    };
    double zeror[10], zeroi[10];
    //rootarr = roots( poly );
    real_poly_roots(poly, 8, zeror, zeroi);
    bigr2= -99999990.0;

    for (int j=1; j <= 8; j++) {
        // if (( rootarr(j) > bigr2 ) && ( isreal(rootarr(j))))  {
        if (( zeror[j] > bigr2 ) && ( zeroi[j]==0))  {
            bigr2= zeror[j];
        }
    }

    u = Constants::GM_Earth/(pow(bigr2,3));
    
    // doubles
    C1 = a1+b1*u;
    C2 = -1;
    C3 = a3+b3*u;

    Matrix aux3(3);
    aux3(1)=C1;
    aux3(2)=C2;
    aux3(3)=C3;

    temp = -D*aux3.transpose();
    rho1 = temp(1)/(a1+b1*u);
    rho2 = -temp(2);
    rho3 = temp(3)/(a3+b3*u);

    rhoold1 = rho1;
    rhoold2 = rho2;
    rhoold3 = rho3;

    rho2 = 99999999.9;
    ll   = 0;

    while ((fabs(rhoold2-rho2) > 1e-12) && (ll <= 0 )) {
        ll = ll + 1;
        rho2 = rhoold2;
        
        r1 = Rs1+Lm1*rho1;
        r2 = Rs2+Lm2*rho2;
        r3 = Rs3+Lm3*rho3;
        
        magr1 = norm(r1);
        magr2 = norm(r2);
        magr3 = norm(r3);
        
        std::tie(v2, theta,theta1,copa,error) = gibbs(r1,r2,r3);
        
        if ( (error!="          ok") & (copa < Constants::pi/180) ) {
            std::tie(v2,theta,theta1,copa,error) = hgibbs(r1,r2,r3,Mjd1,Mjd2,Mjd3);
        }
        Matrix temporal = union_vector(r2.transpose(), v2.transpose());
        auto [p, a, e, i, Omega, omega, M] = elements (temporal);
        
        if ( ll <= 8)  {
            u = Constants::GM_Earth/pow(magr2,3);
            rdot= dot(r2,v2)/magr2;
            udot= (-3*Constants::GM_Earth*rdot)/(pow(magr2,4));
            
            tausqr= tau1*tau1;
            f1=  1 - 0.5*u*tausqr -(1/6)*udot*tausqr*tau1 
                - (1/24) * u*u*tausqr*tausqr 
                - (1/30)*u*udot*tausqr*tausqr*tau1;
            g1= tau1 - (1/6)*u*tau1*tausqr - (1/12) * udot*tausqr*tausqr 
                - (1/120)*u*u*tausqr*tausqr*tau1 
                - (1/120)*u*udot*tausqr*tausqr*tausqr;
            tausqr= tau3*tau3;
            f3=  1 - 0.5*u*tausqr -(1/6)*udot*tausqr*tau3 
                - (1/24) * u*u*tausqr*tausqr 
                - (1/30)*u*udot*tausqr*tausqr*tau3;
            g3= tau3 - (1/6)*u*tau3*tausqr - (1/12) * udot*tausqr*tausqr 
                - (1/120)*u*u*tausqr*tausqr*tau3 
                - (1/120)*u*udot*tausqr*tausqr*tausqr;
        } else {
            
            theta  = angl( r1,r2 );
            theta1 = angl( r2,r3 );
            
            f1= 1 - ( (magr1*(1 - cos(theta)) / p ) );
            g1= ( magr1*magr2*sin(-theta) ) / sqrt( p );
            f3= 1 - ( (magr3*(1 - cos(theta1)) / p ) );
            g3= ( magr3*magr2*sin(theta1) ) / sqrt( p );
        }
        
        C1 = g3/(f1*g3-f3*g1);
        C2 = -1;
        C3 =-g1/(f1*g3-f3*g1);
        
        H1 = Constants::GM_Earth*tau3/12;
        H3 =-Constants::GM_Earth*tau1/12;
        H2 = H1-H3;
        
        G1 = -tau3/(tau1*(tau3-tau1));
        G3 = -tau1/(tau3*(tau3-tau1));
        G2 = G1-G3;
        
        D1 = G1+H1/pow(magr1,3);
        D2 = G2+H2/pow(magr2,3);
        D3 = G3+H3/pow(magr3,3);
        //CX, DX son double
        Matrix aux4(3), aux5(3);
        aux4(1)=D1; aux4(2)=D2; aux4(3)=D3;
        aux5(1)=C1; aux5(2)=C2; aux5(3)=C3;
        double temp2 = (-aux4*aux5.transpose())(1,1);
        rhoold1 = temp2/(a1+b1*u);
        rhoold2 = -temp2;
        rhoold3 = temp2/(a3+b3*u);
        
        r1 = Rs1+Lm1*rhoold1;
        r2 = Rs2+Lm2*rhoold2;
        r3 = Rs3+Lm3*rhoold3;
        
    }

    r1 = Rs1+Lm1*rho1;
    r2 = Rs2+Lm2*rho2;
    r3 = Rs3+Lm3*rho3;

    return std::tie(r2, v2);
}
