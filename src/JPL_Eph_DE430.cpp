#include "..\include\JPL_Eph_DE430.hpp"

Matrix& range(int start, int skip, int end) {
    int num = (end-start)/skip;
    Matrix* m = new Matrix(num+1);
    for(int i = 0; i <=num; i++) {
        (*m)(i+1)= ((double) (start+skip*(i)));
    }

    return *m;
}


std::tuple<Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&, Matrix&> JPL_Eph_DE430(double Mjd_TDB) {
    
    double JD = Mjd_TDB + 2400000.5;
    Matrix temp, PCtemp, Cy_Earth, Cz_Earth, Cx_Earth, Cx, Cy, Cz, Cy_Moon, Cz_Moon, Cx_Moon, Nutations, Cx_Sun, Cz_Sun, Cy_Sun,
    Cx_Mercury, Cy_Mercury, Cz_Mercury, Cx_Venus, Cy_Venus, Cz_Venus, Cx_Mars, Cy_Mars, Cz_Mars, Cx_Jupiter, Cy_Jupiter, Cz_Jupiter, 
    Cx_Saturn, Cy_Saturn, Cz_Saturn, Cx_Uranus, Cy_Uranus, Cz_Uranus, Cx_Neptune,Cy_Neptune,Cz_Neptune, Cx_Pluto,Cy_Pluto,Cz_Pluto,
    Cx_Nutations,Cy_Nutations, Cx_Librations, Cy_Librations, Cz_Librations, Librations;
    
    double i, t1, dt, Mjd0, j, EMRAT, EMRAT1;
    //i = find(PC(:,1)<=JD & JD<=PC(:,2),1,'first');
    for(i=1.0; i < PC.n_row; i++) {
        if (PC(i, 1)<=JD && JD<=PC(i, 2)) break;
    } 
    PCtemp = PC.extract_row(i);//PC(i,:);
    t1 = PCtemp(1)-2400000.5; // MJD at start of interval

    dt = Mjd_TDB - t1;

    //temp = (231:13:270); // 231   244   257   270
    temp = range(231, 13, 270);
    
    Cx_Earth = PCtemp.extract_vector(temp(1), temp(2)-1);
    Cy_Earth = PCtemp.extract_vector(temp(2),temp(3)-1);
    Cz_Earth = PCtemp.extract_vector(temp(3),temp(4)-1);
    
    temp = temp+39;
    Cx = PCtemp.extract_vector(temp(1),temp(2)-1);
    Cy = PCtemp.extract_vector(temp(2),temp(3)-1);
    Cz = PCtemp.extract_vector(temp(3),temp(4)-1);
    Cx_Earth = union_vector(Cx_Earth, Cx);
    Cy_Earth = union_vector(Cy_Earth,Cy);
    Cz_Earth = union_vector(Cz_Earth,Cz);
    if (0<=dt && dt<=16) {
        j=0;
        Mjd0 = t1;
    }
    else if(16<dt && dt<=32) {
        j=1;
        Mjd0 = t1+16*j;
    }
    Matrix &r_Earth = cheb3d(Mjd_TDB, 13, Mjd0, Mjd0+16, Cx_Earth.extract_vector(13*j+1,13*j+13),
                        Cy_Earth.extract_vector(13*j+1,13*j+13), Cz_Earth.extract_vector(13*j+1,13*j+13)).transpose()*1e3;
    
    
    temp = range(441,13,480);
    Cx_Moon = PCtemp.extract_vector(temp(1), temp(2)-1);
    Cy_Moon = PCtemp.extract_vector(temp(2), temp(3)-1);
    Cz_Moon = PCtemp.extract_vector(temp(3), temp(4)-1);
    for (int i = 1.0; i <=7; i++) { //i=1:7
        temp = temp+39;
        Cx = PCtemp.extract_vector(temp(1), temp(2)-1);
        Cy = PCtemp.extract_vector(temp(2), temp(3)-1);
        Cz = PCtemp.extract_vector(temp(3), temp(4)-1);   
        Cx_Moon = union_vector(Cx_Moon, Cx);
        Cy_Moon = union_vector(Cy_Moon, Cy);
        Cz_Moon = union_vector(Cz_Moon, Cz);    
    }
    if (0<=dt && dt<=4) {
        j=0;
        Mjd0 = t1;
    }
    else if(4<dt && dt<=8) {
        j=1;
        Mjd0 = t1+4*j;
    
    } else if(8<dt && dt<=12) {
        j=2;
        Mjd0 = t1+4*j;
    } else if(12<dt && dt<=16) {
        j=3;
        Mjd0 = t1+4*j;
    } else if(16<dt && dt<=20) {
        j=4;
        Mjd0 = t1+4*j;
    } else if(20<dt && dt<=24) {
        j=5;
        Mjd0 = t1+4*j;
    } else if(24<dt && dt<=28) {
        j=6;
        Mjd0 = t1+4*j;
    } else if(28<dt && dt<=32) {
        j=7;
        Mjd0 = t1+4*j;
    }
    Matrix &r_Moon = cheb3d(Mjd_TDB, 13, Mjd0, Mjd0+4, Cx_Moon.extract_vector(13*j+1,13*j+13),
                        Cy_Moon.extract_vector(13*j+1,13*j+13), Cz_Moon.extract_vector(13*j+1,13*j+13)).transpose()*1e3;
    
    temp = range(753, 11, 786);
    Cx_Sun = PCtemp.extract_vector(temp(1), temp(2)-1);
    Cy_Sun = PCtemp.extract_vector(temp(2), temp(3)-1);
    Cz_Sun = PCtemp.extract_vector(temp(3), temp(4)-1);
    temp = temp+33;
    Cx = PCtemp.extract_vector(temp(1), temp(2)-1);
    Cy = PCtemp.extract_vector(temp(2), temp(3)-1);
    Cz = PCtemp.extract_vector(temp(3), temp(4)-1);   
    Cx_Sun = union_vector(Cx_Sun, Cx);
    Cy_Sun = union_vector(Cy_Sun, Cy);
    Cz_Sun = union_vector(Cz_Sun, Cz);
    if (0<=dt && dt<=16) {
        j=0;
        Mjd0 = t1;
    } else if(16<dt && dt<=32) {
        j=1;
        Mjd0 = t1+16*j;
    }
    Matrix &r_Sun = cheb3d(Mjd_TDB, 11, Mjd0, Mjd0+16, Cx_Sun.extract_vector(11*j+1,11*j+11),
                    Cy_Sun.extract_vector(11*j+1,11*j+11), Cz_Sun.extract_vector(11*j+1,11*j+11)).transpose()*1e3;
    temp = range(3, 14, 45);
    Cx_Mercury = PCtemp.extract_vector(temp(1), temp(2)-1);
    Cy_Mercury = PCtemp.extract_vector(temp(2), temp(3)-1);
    Cz_Mercury = PCtemp.extract_vector(temp(3), temp(4)-1);
    for (i=1; i <= 3; i++) {
        temp = temp+42;
        Cx = PCtemp.extract_vector(temp(1), temp(2)-1);
        Cy = PCtemp.extract_vector(temp(2), temp(3)-1);
        Cz = PCtemp.extract_vector(temp(3), temp(4)-1);
        Cx_Mercury = union_vector(Cx_Mercury, Cx);
        Cy_Mercury = union_vector(Cy_Mercury, Cy);
        Cz_Mercury = union_vector(Cz_Mercury, Cz);    
    }
    if (0<=dt && dt<=8) {
        j=0;
        Mjd0 = t1;
    } else if(8<dt && dt<=16) {
        j=1;
        Mjd0 = t1+8*j;
    } else if (16<dt && dt<=24) {
        j=2;
        Mjd0 = t1+8*j;
    } else if(24<dt && dt<=32) {
        j=3;
        Mjd0 = t1+8*j;
    }
    Matrix &r_Mercury = cheb3d(Mjd_TDB, 14, Mjd0, Mjd0+8, Cx_Mercury.extract_vector(14*j+1,14*j+14),
                        Cy_Mercury.extract_vector(14*j+1,14*j+14), Cz_Mercury.extract_vector(14*j+1,14*j+14)).transpose()*1e3;

    temp = range(171, 10, 201);
    Cx_Venus = PCtemp.extract_vector(temp(1), temp(2)-1);
    Cy_Venus = PCtemp.extract_vector(temp(2), temp(3)-1);
    Cz_Venus = PCtemp.extract_vector(temp(3), temp(4)-1);
    temp = temp+30;
    Cx = PCtemp.extract_vector(temp(1), temp(2)-1);
    Cy = PCtemp.extract_vector(temp(2), temp(3)-1);
    Cz = PCtemp.extract_vector(temp(3), temp(4)-1);
    Cx_Venus = union_vector(Cx_Venus, Cx);
    Cy_Venus = union_vector(Cy_Venus, Cy);
    Cz_Venus = union_vector(Cz_Venus, Cz);
    if (0<=dt && dt<=16) {
        j=0;
        Mjd0 = t1;
    } else if(16<dt && dt<=32) {
        j=1;
        Mjd0 = t1+16*j;
    }
    Matrix &r_Venus = cheb3d(Mjd_TDB, 10, Mjd0, Mjd0+16, Cx_Venus.extract_vector(10*j+1,10*j+10),
                        Cy_Venus.extract_vector(10*j+1,10*j+10), Cz_Venus.extract_vector(10*j+1,10*j+10)).transpose()*1e3;

    temp = range(309, 11, 342);
    Cx_Mars = PCtemp.extract_vector(temp(1), temp(2)-1);
    Cy_Mars = PCtemp.extract_vector(temp(2), temp(3)-1);
    Cz_Mars = PCtemp.extract_vector(temp(3), temp(4)-1);
    j=0;
    Mjd0 = t1;
    Matrix &r_Mars = cheb3d(Mjd_TDB, 11, Mjd0, Mjd0+32, Cx_Mars.extract_vector(11*j+1,11*j+11),
                        Cy_Mars.extract_vector(11*j+1,11*j+11), Cz_Mars.extract_vector(11*j+1,11*j+11)).transpose()*1e3;

    temp = range(342, 8, 366);
    Cx_Jupiter = PCtemp.extract_vector(temp(1), temp(2)-1);
    Cy_Jupiter = PCtemp.extract_vector(temp(2), temp(3)-1);
    Cz_Jupiter = PCtemp.extract_vector(temp(3), temp(4)-1);
    j=0;
    Mjd0 = t1;
    Matrix &r_Jupiter = cheb3d(Mjd_TDB, 8, Mjd0, Mjd0+32, Cx_Jupiter.extract_vector(8*j+1,8*j+8),
                        Cy_Jupiter.extract_vector(8*j+1,8*j+8), Cz_Jupiter.extract_vector(8*j+1,8*j+8)).transpose()*1e3;

    temp = range(366, 7, 387);
    Cx_Saturn = PCtemp.extract_vector(temp(1), temp(2)-1);
    Cy_Saturn = PCtemp.extract_vector(temp(2), temp(3)-1);
    Cz_Saturn = PCtemp.extract_vector(temp(3), temp(4)-1);
    j=0;
    Mjd0 = t1;
    Matrix &r_Saturn = cheb3d(Mjd_TDB, 7, Mjd0, Mjd0+32, Cx_Saturn.extract_vector(7*j+1,7*j+7),
                        Cy_Saturn.extract_vector(7*j+1,7*j+7), Cz_Saturn.extract_vector(7*j+1,7*j+7)).transpose()*1e3;

    temp = range(387, 6, 405);
    Cx_Uranus = PCtemp.extract_vector(temp(1), temp(2)-1);
    Cy_Uranus = PCtemp.extract_vector(temp(2), temp(3)-1);
    Cz_Uranus = PCtemp.extract_vector(temp(3), temp(4)-1);
    j=0;
    Mjd0 = t1;
    Matrix &r_Uranus = cheb3d(Mjd_TDB, 6, Mjd0, Mjd0+32, Cx_Uranus.extract_vector(6*j+1,6*j+6),
                        Cy_Uranus.extract_vector(6*j+1,6*j+6), Cz_Uranus.extract_vector(6*j+1,6*j+6)).transpose()*1e3;

    temp = range(405, 6, 423);
    Cx_Neptune = PCtemp.extract_vector(temp(1), temp(2)-1);
    Cy_Neptune = PCtemp.extract_vector(temp(2), temp(3)-1);
    Cz_Neptune = PCtemp.extract_vector(temp(3), temp(4)-1);
    j=0;
    Mjd0 = t1;
    Matrix &r_Neptune = cheb3d(Mjd_TDB, 6, Mjd0, Mjd0+32, Cx_Neptune.extract_vector(6*j+1,6*j+6),
                        Cy_Neptune.extract_vector(6*j+1,6*j+6), Cz_Neptune.extract_vector(6*j+1,6*j+6)).transpose()*1e3;

    temp = range(423, 6, 441);
    Cx_Pluto = PCtemp.extract_vector(temp(1), temp(2)-1);
    Cy_Pluto = PCtemp.extract_vector(temp(2), temp(3)-1);
    Cz_Pluto = PCtemp.extract_vector(temp(3), temp(4)-1);
    j=0;
    Mjd0 = t1;
    Matrix &r_Pluto = cheb3d(Mjd_TDB, 6, Mjd0, Mjd0+32, Cx_Pluto.extract_vector(6*j+1,6*j+6),
                        Cy_Pluto.extract_vector(6*j+1,6*j+6), Cz_Pluto.extract_vector(6*j+1,6*j+6)).transpose()*1e3;

    temp = range(819, 10, 839);
    Cx_Nutations = PCtemp.extract_vector(temp(1), temp(2)-1);
    Cy_Nutations = PCtemp.extract_vector(temp(2), temp(3)-1);
    for (i=1; i <= 3; i++) {
        temp = temp+20;
        Cx = PCtemp.extract_vector(temp(1), temp(2)-1);
        Cy = PCtemp.extract_vector(temp(2), temp(3)-1);
        Cx_Nutations = union_vector(Cx_Nutations, Cx);
        Cy_Nutations = union_vector(Cy_Nutations, Cy);
    }
    if (0<=dt && dt<=8) {
        j=0;
        Mjd0 = t1;
    } else if(8<dt && dt<=16) {
        j=1;
        Mjd0 = t1+8*j;
    } else if (16<dt && dt<=24) {
        j=2;
        Mjd0 = t1+8*j;
    } else if(24<dt && dt<=32) {
        j=3;
        Mjd0 = t1+8*j;
    }
    Nutations = cheb3d(Mjd_TDB, 10, Mjd0, Mjd0+8, Cx_Nutations.extract_vector(10*j+1,10*j+10),
                    Cy_Nutations.extract_vector(10*j+1,10*j+10),zeros(10,1)).transpose();

    temp = range(899, 10, 929);
    Cx_Librations = PCtemp.extract_vector(temp(1), temp(2)-1);
    Cy_Librations = PCtemp.extract_vector(temp(2), temp(3)-1);
    Cz_Librations = PCtemp.extract_vector(temp(3), temp(4)-1);
    for (i=1; i <= 3; i++) {
        temp = temp+30;
        Cx = PCtemp.extract_vector(temp(1), temp(2)-1);
        Cy = PCtemp.extract_vector(temp(2), temp(3)-1);
        Cz = PCtemp.extract_vector(temp(3), temp(4)-1);
        Cx_Librations = union_vector(Cx_Librations, Cx);
        Cy_Librations = union_vector(Cy_Librations, Cy);
        Cz_Librations = union_vector(Cz_Librations, Cz);    
    }
    if (0<=dt && dt<=8) {
        j=0;
        Mjd0 = t1;
    } else if(8<dt && dt<=16) {
        j=1;
        Mjd0 = t1+8*j;
    } else if (16<dt && dt<=24) {
        j=2;
        Mjd0 = t1+8*j;
    } else if(24<dt && dt<=32) {
        j=3;
        Mjd0 = t1+8*j;
    }
    Librations = cheb3d(Mjd_TDB, 10, Mjd0, Mjd0+8, Cx_Librations.extract_vector(10*j+1,10*j+10),
                        Cy_Librations.extract_vector(10*j+1,10*j+10), Cz_Librations.extract_vector(10*j+1,10*j+10)).transpose();
    EMRAT = 81.30056907419062; // DE430
    EMRAT1 = 1/(1+EMRAT);
    r_Earth = r_Earth-r_Moon*EMRAT1;
    r_Mercury = -r_Earth+r_Mercury;
    r_Venus = -r_Earth+r_Venus;
    r_Mars = -r_Earth+r_Mars;
    r_Jupiter = -r_Earth+r_Jupiter;
    r_Saturn = -r_Earth+r_Saturn;
    r_Uranus = -r_Earth+r_Uranus;
    r_Neptune = -r_Earth+r_Neptune;
    r_Pluto = -r_Earth+r_Pluto;
    r_Sun = -r_Earth+r_Sun;

    return std::tie(r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, r_Neptune,r_Pluto,r_Moon,r_Sun);


}
