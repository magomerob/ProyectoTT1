// $Source$
/**
 * @file accelHarmonic.cpp
 * @author Marcos Gómez Robres
 * @brief Implementación de la aceleración armónica.
 * @version 0.1
 * @date 2025-04-30
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */


#include "../include/accelHarmonic.h"

 Matrix accelHarmonic (Matrix &r, Matrix &E, double n_max, double m_max){
    Matrix Cnm = zeros(181,181);
    Matrix Snm = zeros(181,181);
    FILE *fid = fopen("../data/GGM03S.txt", "r");
    double aux;
    for(int i = 1; i <= 181; i++){
        for(int j = 1; j <= i; j++){
            fscanf(fid, "%lf %lf %lf %lf %lf %lf",&aux, &aux, &(Cnm(i,j)), &(Snm(i,j)), &aux, &aux);
        }
    }
    
    double r_ref = 6378.1363e3;
    double gm    = 398600.4415e9;

    Matrix r_bf = E * r;
    double d = norm(r_bf);
    //cout<<"d"<<d<<endl;
    double latgc = asin(r_bf(3,1)/d*1.0);
    //cout<<"latgc"<<latgc<<endl;
    double lon = atan2(r_bf(2,1)*1.0,r_bf(1,1)*1.0);
    //cout<<"lon"<<lon<<endl;
    auto [pnm, dpnm] = legendre(n_max,m_max,latgc);
    
    double dUdr = 0;
    double dUdlatgc = 0;
    double dUdlon = 0;
    double q3 = 0;
    double q2 = 0;
    double q1 = 0;
    
    for (int n = 0; n <= n_max; n++) {
        double b1 = ((-gm) / (d * d)) * pow((r_ref / d), n) * (n + 1); 
        double b2 = (gm / d) * pow((r_ref / d), n);
        double b3 = (gm / d) * pow((r_ref / d), n);
        
        for (int m = 0; m <= m_max; m++) {
            q1 += pnm(n + 1, m + 1) * (Cnm(n + 1, m + 1) * cos(1.0*m * lon) + Snm(n + 1, m + 1) * sin(1.0*m * lon));
            q2 += dpnm(n + 1, m + 1) * (Cnm(n + 1, m + 1) * cos(1.0*m * lon) + Snm(n + 1, m + 1) * sin(1.0*m * lon));
            q3 += m * pnm(n + 1, m + 1) * (Snm(n + 1, m + 1) * cos(1.0*m * lon) - Cnm(n + 1, m + 1) * sin(1.0*m * lon));
        }
    
        dUdr     += q1 * b1;
        dUdlatgc += q2 * b2;
        dUdlon   += q3 * b3;
    
        q1 = 0; q2 = 0; q3 = 0;

        

    }
    
    double r2xy = pow(r_bf(1, 1), 2) + pow(r_bf(2, 1), 2);
    double sqrt_r2xy = sqrt(r2xy);
    
    double ax = (1 / d * dUdr - r_bf(3, 1) / (d * d * sqrt_r2xy) * dUdlatgc) * r_bf(1, 1)
               - (1 / r2xy * dUdlon) * r_bf(2, 1);
    
    double ay = (1 / d * dUdr - r_bf(3, 1) / (d * d * sqrt_r2xy) * dUdlatgc) * r_bf(2, 1)
               + (1 / r2xy * dUdlon) * r_bf(1, 1);
    
    double az = (1 / d * dUdr) * r_bf(3, 1) + (sqrt_r2xy / (d * d)) * dUdlatgc;

    Matrix a_bf(3);
    a_bf(1) = ax;
    a_bf(2) = ay;
    a_bf(3) = az;

    Matrix a = transpose(E) * transpose(a_bf);
    return a;
 }