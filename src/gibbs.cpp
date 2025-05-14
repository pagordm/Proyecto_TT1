#include "..\include\gibbs.hpp"

/**
 * @brief Performs the Gibbs method of orbit determination to find velocity at middle point of 3 position vectors
 * 
 * This method determines the velocity at the middle point of the 3 given position vectors.
 * 
 * @param r1 IJK position vector #1 (m)
 * @param r2 IJK position vector #2 (m) 
 * @param r3 IJK position vector #3 (m)
 * @return std::tuple<Matrix&,double,std::string> Returns:
 *         - v2: IJK velocity vector for r2 (m/s)
 *         - theta: Angle between vectors (rad)
 *         - error: Flag indicating success ('ok',...)
 */
std::tuple<Matrix&, double, double, double, std::string> gibbs(Matrix r1, Matrix r2, Matrix r3) {
    double small, theta, theta1, magr1, magr2, magr3, copa, magd, magn, r1mr2, r3mr1, r2mr3, l, tover2;
    Matrix *v2 = new Matrix(3); Matrix p, q, w, pn, r1n, d, n, nn, dn, s, b;
    std::string error;
    
    small=0.00000001;
    theta= 0.0;
    error = "          ok";
    theta1= 0.0;
    magr1 = norm( r1.transpose() );
    magr2 = norm( r2.transpose() );
    magr3 = norm( r3.transpose() );
    for (int i= 1; i <=3; i++) {
        (*v2)(i)= 0.0;
    }

    p = cross( r2.transpose(),r3.transpose() );
    q = cross( r3.transpose(),r1.transpose() );
    w = cross( r1.transpose(),r2.transpose() );
    pn = unit( p.transpose() );
    r1n = unit( r1 );
    copa=  asin( dot( pn,r1n ) );

    if ( abs( dot(r1n.transpose(),pn.transpose()) ) > 0.017452406 )  {
        error= "not coplanar";
    }

    d = p + q + w;
    magd = norm(d);
    n = p*magr1 + q*magr2 + w*magr3;
    magn = norm(n);
    nn = unit( n );
    dn = unit( d );

    // -------------------------------------------------------------
    // determine if  the orbit is possible. both d and n must be in
    // the same direction, and non-zero.
    // -------------------------------------------------------------
    if ( ( abs(magd)<small ) || ( abs(magn)<small ) || ( dot(nn,dn) < small ) ) {
        error= "  impossible";
    } else {
        theta  = angl( r1.transpose(),r2.transpose() );
        //cout << "gibbs theta: " << theta << endl;
        theta1 = angl( r2.transpose(),r3.transpose() );
        //cout << "gibbs theta1: " << theta1 << endl;

        // ----------- perform gibbs method to find v2 -----------
        r1mr2= magr1-magr2;
        r3mr1= magr3-magr1;
        r2mr3= magr2-magr3;
        s  = r3*r1mr2 + r2*r3mr1 + r1*r2mr3;
        b  = cross( d,r2.transpose() );
        l  = sqrt(Constants::GM_Earth / (magd*magn) );
        tover2= l / magr2;
        *v2 = b.transpose()*tover2 + s*l;
    }

    return std::tie(*v2, theta, theta1, copa, error);
}