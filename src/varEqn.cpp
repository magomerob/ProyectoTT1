// $Source$
/**
 * @file varEqn.cpp
 * @author Marcos Gómez Robres
 * @brief Implementación de varEqn.
 * @version 0.1
 * @date 2025-05-14
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#include "../include/varEqn.h"

Matrix varEqn(double x, Matrix yPhi)
{

    Matrix tyPhi  = transpose(yPhi);

    auto [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = iers(eopdata,AuxParam.Mjd_UTC,'l');
    auto [UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
    double Mjd_UT1 = AuxParam.Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;

    // Transformation matrix
    Matrix P = precMatrix(SAT_Const::MJD_J2000,AuxParam.Mjd_TT + x/86400.0);
    Matrix N = nutMatrix(AuxParam.Mjd_TT + x/86400.0);
    Matrix T = N * P;
    Matrix E = poleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    // State vector components
    Matrix r = extract_vector(tyPhi,1,3);
    Matrix v = extract_vector(tyPhi,4,6);
    Matrix Phi = zeros(6,6);

    // State transition matrix
    for(size_t j=1; j<=6; j++){
        Phi.assign_column(extract_vector(tyPhi,6*j+1,6*j+6),j);
    }
    
    // Acceleration and gradient
    Matrix a = accelHarmonic ( transpose(r), E, AuxParam.n, AuxParam.m );
    Matrix G = g_accelHarmonic ( transpose(r), E, AuxParam.n, AuxParam.m );

    // Time derivative of state transition matrix
    Matrix yPhip = zeros(42,1);
    Matrix dfdy = zeros(6,6);

    for (size_t i=1; i<=3; i++) {
        for(size_t j=1; j<=3; j++) {
            dfdy(i,j) = 0.0;                 // dv/dr(i,j)
            dfdy(i+3,j) = G(i,j);            // da/dr(i,j)
            if ( i==j ){
                dfdy(i,j+3) = 1;
            }else{
                dfdy(i,j+3) = 0;             // dv/dv(i,j)
            }
            dfdy(i+3,j+3) = 0.0;             // da/dv(i,j)
        }
    }


    Matrix Phip = dfdy*Phi;

    // Derivative of combined state vector and state transition matrix
    for (size_t i=1; i<=3; i++) {
        yPhip(i,1)   = v(i);                 // dr/dt(i)
        yPhip(i+3,1) = a(i,1);                 // dv/dt(i)
    }

    for (size_t j=1; j<=6; j++) {
        for (size_t i=1; i<=6; i++) {
            yPhip(6*j+i,1) = Phip(i,j);     // dPhi/dt(i,j)
        }
    }

    return yPhip;
}