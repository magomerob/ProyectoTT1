// $Source$
/**
 * @file global.cpp
 * @author Marcos Gómez Robres
 * @brief Esta es la implementación de la clase global.
 * @date 2025-04-23
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

 
 #include "../include/global.h"
 

Matrix eopdata;

void eop19620101(int c){
    eopdata = Matrix(c,13);
    FILE *fid = fopen("../data/eop19620101.txt", "r");
    if(fid == NULL){
        printf("Fail opening eop19620101.txt\n");
        exit(EXIT_FAILURE);
    }
    for(int i = 1; i<=c; i++){

        fscanf(fid,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&eopdata(i,1),&eopdata(i,2),&eopdata(i,3),&eopdata(i,4),&eopdata(i,5),&eopdata(i,6),&eopdata(i,7),&eopdata(i,8),&eopdata(i,9),&eopdata(i,10),&eopdata(i,11),&eopdata(i,12),&eopdata(i,13));
    }
	/*for(int i = 1; i<=13; i++){
		cout<<eopdata(1,i)<<" ";
	}*/
    fclose(fid);
}

Matrix PC;

void DE430Coeff() {
    PC = Matrix(2285, 1020);
	FILE *fid = fopen("../data/DE430Coeff.txt","r");

	if(fid== NULL) {
		printf("Can't DE430Coeff.txt\n");
		exit(EXIT_FAILURE);
	}
	double aux;
	for (int n = 1; n <= 2285; n++) {
		for(int m=1;m<=1020;m++){
				fscanf(fid, "%lf,",&(PC(n, m)));
			}
		}
	fclose(fid);
}

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