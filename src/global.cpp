#include "..\include\global.hpp"
#include "..\include\matrix.h"


Matrix eopdata;
Param AuxParam;

void auxparam() {
    AuxParam.Mjd_UTC = 49746.1163541665;
    AuxParam.n      = 20;
    AuxParam.m      = 20;
    AuxParam.sun     = 1;
    AuxParam.moon    = 1;
    AuxParam.planets = 1;
    AuxParam.Mjd_TT  = 49746.1170623147;
}

void eop19620101(int c) {
    eopdata = zeros(13, c);

    FILE *fid = fopen("..\\data\\eop19620101.txt", "r");
    if (fid== NULL) {
        printf("Fail open eop19620101.txt file\n");
        exit(EXIT_FAILURE);
    }
    for(int i = 1; i <= c; i++) {
        fscanf(fid,"%lf %lf %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            &(eopdata(1,i)),&(eopdata(2,i)),&(eopdata(3,i)),&(eopdata(4,i)),&(eopdata(5,i)),
            &(eopdata(6,i)),&(eopdata(7,i)),&(eopdata(8,i)),&(eopdata(9,i)),&(eopdata(10,i)),
            &(eopdata(11,i)),&(eopdata(12,i)),&(eopdata(13,i)));
    }
    fclose(fid);
}

Matrix PC;

void DE430Coeff() {
    PC = zeros(2285, 1020);
	FILE *fid = fopen("../data/DE430Coeff.txt","r");

	if(fid== NULL) {
		printf("Fail open DE430Coeff.txt file\n");
		exit(EXIT_FAILURE);
	}
	double aux;
	for (int n = 1; n <= 2285; n++) {
		for(int m=1;m<=1020;m++){
				fscanf(fid, "%lf ",&(PC(n, m)));
			}
		}
	fclose(fid);
}


Matrix Cnm;
Matrix Snm;

void GGM03S() {
    
    Cnm = zeros(181, 181);
    Snm = zeros(181, 181);
    
    FILE *fid = fopen("..\\data\\GGM03S.txt", "r");
    if (fid== NULL) {
        printf("Fail open GGM03S.txt file\n");
        exit(EXIT_FAILURE);
    }

    double aux;
    for(int n = 0; n <= 180; n++) {
        for(int m = 0; m <= n; m++) {
            fscanf(fid,"%lf %lf %lf %lf %lf %lf",
                &aux,&aux,&(Cnm(n+1, m+1)),&(Snm(n+1, m+1)),
                &aux,&aux);
        }
    }
    fclose(fid);
}

Matrix obs;

void GEOS3(int f) {
    obs = zeros(f, 4);

    FILE *fp = fopen("../data/GEOS3.txt", "r");

    int Y, MO, D, H, MI, S;
    double AZ, EL, DIST;
    char line[55], y[5], mo[3], d[3], h[3], mi[3], s[3], az[9], el[9], dist[10];
    for(int i = 1; i <=f; i++) {
        fgets(line, sizeof(line)+2, fp);

        strncpy(y, &(line[0]), 4);
        y[4]='\0';
        Y=atoi(y);

        strncpy(mo, &(line[5]), 2);
        mo[2]='\0';
        MO = atoi(mo);

        strncpy(d, &(line[8]), 2);
        d[2]='\0';
        D = atoi(d);

        strncpy(h, &(line[12]), 2);
        h[2]='\0';
        H = atoi(h);

        strncpy(mi, &(line[15]), 2);
        mi[2]='\0';
        MI=atoi(mi);

        strncpy(s, &(line[18]), 2);
        s[2]='\0';
        S = atoi(s);

        strncpy(az, &(line[25]), 8);
        az[8]='\0';
        AZ = atof(az);

        strncpy(el, &(line[34]), 8);
        el[8]='\0';
        EL=atof(el);

        strncpy(dist, &(line[43]), 9);
        dist[9]='\0';
        DIST = atof(dist);

        obs(i, 1) = Mjday(Y, MO, D, H, MI, S);
        obs(i, 2) = Constants::Rad*AZ;
        obs(i, 3) = Constants::Rad*EL;
        obs(i, 4) = 1e3*DIST;

    }
}
