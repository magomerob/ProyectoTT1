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

    //cout<<"r:"<<endl<<r<<endl;
    //cout<<"E:"<<endl<<E<<endl;

    double r_ref = 6378.1363e3;
    double gm    = 398600.4415e9;

    Matrix r_bf = E * r;

    //cout<<"r_bf:"<<endl<<r_bf<<endl;

    double d = norm(r_bf)*1.0;

    double latgc = asin(r_bf(3,1)/d);
    //cout<<"latgc"<<latgc<<endl;
    double lon = atan2(r_bf(2,1),r_bf(1,1));
    //cout<<"lon"<<lon<<endl;
    auto [pnm, dpnm] = legendre(n_max,m_max,latgc);

    //cout<<"pnm: "<<pnm<<endl;
    //cout<<"dpnm: "<<dpnm<<endl;

    long double dUdr = 0;
    long double dUdlatgc = 0;
    long double dUdlon = 0;
    double q3 = 0.0;
    double q2 = 0.0;
    double q1 = 0.0;
    
    double b1,b2,b3;

    for (int n = 0; n <= n_max; n++) {
        b1 = (-gm/pow(d, 2)*1.0)*pow((r_ref/d),n)*(n+1);
        b2 = (gm / d) * pow((r_ref / d), n);
        b3 = (gm / d) * pow((r_ref / d), n);
        
        for (double m = 0; m <= m_max; m++) {
            q1 = q1 + pnm(n+1,m+1)*(Cnm(n+1,m+1)*cos(m*lon)+Snm(n+1,m+1)*sin(m*lon));
            q2 = q2 + dpnm(n+1,m+1)*(Cnm(n+1,m+1)*cos(m*lon)+Snm(n+1,m+1)*sin(m*lon));
            q3 = q3 + m*pnm(n+1,m+1)*(Snm(n+1,m+1)*cos(m*lon)-Cnm(n+1,m+1)*sin(m*lon));
        }
    
        dUdr     = dUdr + q1 * b1;
        dUdlatgc = dUdlatgc + q2 * b2;
        dUdlon   = dUdlon + q3 * b3;
    
        q1 = 0; q2 = 0; q3 = 0;

    }
    
    //cout<<"dUdr: "<<dUdr<<endl;
    //cout<<"dUdlatgc: "<<dUdlatgc<<endl;
    //cout<<"dUdlon: "<<dUdlon<<endl;

    double r2xy = pow(r_bf(1, 1), 2) + pow(r_bf(2, 1), 2);

    double ax = (1.0/d*dUdr-r_bf(3)/(pow(d,2)*sqrt(r2xy))*dUdlatgc)*r_bf(1)-(1.0/r2xy*dUdlon)*r_bf(2);
    double ay = (1.0/d*dUdr-r_bf(3)/(pow(d,2)*sqrt(r2xy))*dUdlatgc)*r_bf(2)+(1.0/r2xy*dUdlon)*r_bf(1);
    double az =  1.0/d*dUdr*r_bf(3)+sqrt(r2xy)/pow(d,2)*dUdlatgc;

    Matrix a_bf(3);
    a_bf(1) = ax;
    a_bf(2) = ay;
    a_bf(3) = az;

    Matrix a = transpose(E) * transpose(a_bf);
    return a;
 }