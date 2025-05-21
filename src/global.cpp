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
		printf("Can't open DE430Coeff.txt\n");
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
    AuxParam.Mjd_UTC = 4.974611635416653e+04;
    AuxParam.n      = 20;
    AuxParam.m      = 20;
    AuxParam.sun     = 1;
    AuxParam.moon    = 1;
    AuxParam.planets = 1;
    AuxParam.Mjd_TT  = 4.974611706231468e+04;
}

Matrix Cnm;
Matrix Snm;

void GGM03S()
{
    Cnm = Matrix(181,181);
    Snm = Matrix(181,181);
    FILE *fid = fopen("../data/GGM03S.txt","r");
    if (fid== NULL) {
        printf("Can't open GGM03S.txt\n");
        exit(EXIT_FAILURE);
    }
    double _;
    for(int n=0; n<=180; n++){
        for(int m = 0; m<=n; m++){
            fscanf(fid,"%lf %lf %lf %lf %lf %lf",&_,&_,&(Cnm(n+1, m+1)),&(Snm(n+1, m+1)),&_,&_);
        }
    }
    fclose(fid);
}

Matrix obs;

void GEOS3()
{
    int nobs = 46;
    obs = Matrix(nobs,4);

    // read observations
    FILE *fid = fopen("../data/GEOS3.txt","r");
    
    int y, mo, d, h, mi, s;
    double az, el, dist;
    char tline[55], y_[5], mo_[3], d_[3], h_[3], mi_[3], s_[3], az_[9], el_[9], dist_[10];

    for(int i= 1; i<=nobs; i++){
        fgets(tline, sizeof(tline)+2, fid);

        strncpy(y_, &(tline[0]), 4);
        y_[4]='\0';
        y=atoi(y_);
        strncpy(mo_, &(tline[5]), 2);
        mo_[2]='\0';
        mo = atoi(mo_);
        strncpy(d_, &(tline[8]), 2);
        d_[2]='\0';
        d = atoi(d_);
        strncpy(h_, &(tline[12]), 2);
        h_[2]='\0';
        h = atoi(h_);
        strncpy(mi_, &(tline[15]), 2);
        mi_[2]='\0';
        mi=atoi(mi_);
        strncpy(s_, &(tline[18]), 2);
        s_[2]='\0';
        s = atoi(s_);
        strncpy(az_, &(tline[25]), 8);
        az_[8]='\0';
        az = atof(az_);
        strncpy(el_, &(tline[34]), 8);
        el_[8]='\0';
        el=atof(el_);
        strncpy(dist_, &(tline[43]), 9);
        dist_[9]='\0';
        dist = atof(dist_);

        obs(i, 1) = Mjday(y, mo, d, h, mi, s);
        obs(i, 2) = SAT_Const::Rad*az;
        obs(i, 3) = SAT_Const::Rad*el;
        obs(i, 4) = 1e3*dist;
    }

    fclose(fid);
}