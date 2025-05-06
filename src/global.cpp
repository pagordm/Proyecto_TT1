#include "..\include\global.hpp"
#include "..\include\matrix.h"

Matrix eopdata;

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