/**
 * @file EKF_GEOS3.cpp
 * @author Marcos Gómez Robres
 * @brief Archivo principal
 * @version 0.1
 * @date 2025-05-18
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#include "../include/matrix.h"
#include "../include/global.h"
#include "../include/position.h"
#include "../include/mjday.h"
#include "../include/accel.h"
#include "../include/LTC.h"
#include "../include/iers.h"
#include "../include/timediff.h"
#include "../include/varEqn.h"
#include "../include/r_z.h"
#include "../include/azelpa.h"
#include "../include/measUpdate.h"
#include "../include/DEInteg.h"
#include "..\include\timeUpdate.h"
#include <iostream>

using namespace std;

int main(){
    cout<<"Running EKF_GEOS3..."<<endl;
    eop19620101(21413);
    DE430Coeff();
    auxparam();
    GGM03S();
    GEOS3();

    double sigma_range, sigma_az, sigma_el,lat,lon,alt,Mjd1,Mjd2,Mjd3,Mjd0,Mjd_UTC,Mjd_TT,Mjd_UT1,theta, Dist;
    double n,m,n_eqn,t, t_old;
    Matrix Y_old,K;

    sigma_range = 92.5;          // [m]
    sigma_az = 0.0224*SAT_Const::Rad; // [rad]
    sigma_el = 0.0139*SAT_Const::Rad; // [rad]

    // Kaena Point station
    lat = SAT_Const::Rad*21.5748;     // [rad]
    lon = SAT_Const::Rad*(-158.2706); // [rad]
    alt = 300.20;                // [m]

    Matrix Rs = transpose(Position(lon, lat, alt));

    Mjd1 = obs(1,1);
    Mjd2 = obs(9,1);
    Mjd3 = obs(18,1);

    Matrix r2(3,1);
    r2(1,1) = 6.221397628578685e+06;
    r2(2,1) = 2.867713779657379e+06;
    r2(3,1) = 3.006155985099489e+06;
    Matrix v2(3,1);
    v2(1,1) = 4.645047251618060e+03;
    v2(2,1) = -2.752215915882042e+03;
    v2(3,1) = -7.507999409870306e+03;

    /*
    [r2,v2] = anglesg(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),...
                    Mjd1,Mjd2,Mjd3,Rs,Rs,Rs);
    % [r2,v2] = anglesdr(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),...
    %                    Mjd1,Mjd2,Mjd3,Rs,Rs,Rs);
    */
    Matrix Y0_apr = transpose(union_vector(transpose(r2),transpose(v2)));

    Mjd0 = Mjday(1995,1,29,02,38,0);

    Mjd_UTC = obs(9,1);

    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.n      = 20;
    AuxParam.m      = 20;
    AuxParam.sun     = 1;
    AuxParam.moon    = 1;
    AuxParam.planets = 1;

    n_eqn  = 6;

    Matrix Y = DEInteg(accel,0,-(obs(9,1)-Mjd0)*86400.0,1e-13,1e-6,6,Y0_apr);

    Matrix P(6,6);
    
    for(int i=1; i<=3; i++){
        P(i,i)=1e8;
    }
    for(int i=4; i<=6; i++){
        P(i,i)=1e3;
    }

    Matrix LT = LTC(lon,lat);

    Matrix yPhi(42,1);
    Matrix Phi(6,6);

    // Measurement loop
    t = 0;
    int nobs = 46;

    for(int i=1; i<=nobs; i++){ 
        // Previous step
        t_old = t;
        Y_old = Y;
        
        // Time increment and propagation
        Mjd_UTC = obs(i,1);                       // Modified Julian Date
        t       = (Mjd_UTC-Mjd0)*86400.0;         // Time since epoch [s]
        
        auto [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = iers(eopdata,Mjd_UTC,'l');
        auto [UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
        Mjd_TT = Mjd_UTC + TT_UTC/86400;
        Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;
        AuxParam.Mjd_UTC = Mjd_UTC;
        AuxParam.Mjd_TT = Mjd_TT;
            
        for(int ii=1; ii<=6; ii++){
            yPhi(ii,1) = Y_old(ii,1);
            for(int j=1; j<=6; j++){ 
                if (ii==j){
                    yPhi(6*j+ii,1) = 1; 
                }else{
                    yPhi(6*j+ii,1) = 0;
                }
            }
        }
        
        yPhi = DEInteg (varEqn,0,t-t_old,1e-13,1e-6,42,yPhi);
        
        // Extract state transition matrices
        for(int j=1; j<=6; j++){
            Phi.assign_column(extract_vector(transpose(yPhi),6*j+1,6*j+6),j);
        }
        
        Y = DEInteg (accel,0,t-t_old,1e-13,1e-6,6,Y_old);
        
        // Topocentric coordinates
        theta = gmst(Mjd_UT1);                    // Earth rotation
        Matrix U = r_z(theta);
        Matrix r = transpose(extract_vector(transpose(Y),1,3));
        Matrix s = LT*(U*r-Rs);                          // Topocentric position [m]
        
        // Time update
        P = TimeUpdate(P, Phi);
            
        // Azimuth and partials
        auto [Azim, Elev, dAds, dEds] = azelpa(s);     // Azimuth, Elevation
        Matrix tmp = dAds*LT*U;
        Matrix dAdY = union_vector(tmp,zeros(1,3));
        
        // Measurement update
        tie(K, Y, P) = measUpdate ( Y, obs(i,2), Azim, sigma_az, dAdY, P, 6 );
        
        // Elevation and partials
        r = transpose(extract_vector(transpose(Y),1,3));
        s = LT*(U*r-Rs);                          // Topocentric position [m]
        tie(Azim, Elev, dAds, dEds) = azelpa(s);     // Azimuth, Elevation
        tmp = dEds*LT*U;
        Matrix dEdY = union_vector(tmp,zeros(1,3));
        
        // Measurement update
        tie(K, Y, P) = measUpdate ( Y, obs(i,3), Elev, sigma_el, dEdY, P, 6 );
        
        // Range and partials
        r = transpose(extract_vector(transpose(Y),1,3));
        s = LT*(U*r-Rs);                          // Topocentric position [m]
        Dist = norm(s);
        Matrix dDds = transpose(s/Dist);         // Range
        tmp = dDds*LT*U;
        Matrix dDdY = union_vector(tmp,zeros(1,3));
        
        // Measurement update
        tie(K, Y, P) = measUpdate ( Y, obs(i,4), Dist, sigma_range, dDdY, P, 6 );
    }

    auto [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = iers(eopdata,obs(46,1),'l');
    auto [UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC/86400;
    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.Mjd_TT = Mjd_TT;

    Matrix Y0 = DEInteg (accel,0,-(obs(46,1)-obs(1,1))*86400.0,1e-13,1e-6,6,Y);

    Matrix Y_true(6,1);
    Y_true(1,1) = 5753.173e3;
    Y_true(2,1) = 2673.361e3;
    Y_true(3,1) = 3440.304e3;
    Y_true(4,1) = 4.324207e3;
    Y_true(5,1) = -1.924299e3;
    Y_true(6,1) = -5.728216e3;

    printf("\nError of Position Estimation\n");
    printf("dX%10.1f [m]\n",Y0(1)-Y_true(1));
    printf("dY%10.1f [m]\n",Y0(2)-Y_true(2));
    printf("dZ%10.1f [m]\n",Y0(3)-Y_true(3));
    printf("\nError of Velocity Estimation\n");
    printf("dVx%8.1f [m/s]\n",Y0(4)-Y_true(4));
    printf("dVy%8.1f [m/s]\n",Y0(5)-Y_true(5));
    printf("dVz%8.1f [m/s]\n",Y0(6)-Y_true(6));
}
