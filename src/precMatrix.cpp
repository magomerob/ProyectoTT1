// $Source$
/**
 * @file precMatrix.cpp
 * @author Marcos Gómez Robres
 * @brief Implementación de la transformación de la matriz de precesión.
 * @version 0.1
 * @date 2025-05-07
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

 #include "../include/precMatrix.h"

 Matrix precMatrix(double Mjd_1, double Mjd_2)
 {
    double T  = (Mjd_1-SAT_Const::MJD_J2000)/36525;
    double dT = (Mjd_2-Mjd_1)/36525;

    // Precession angles
    double zeta  =  ( (2306.2181+(1.39656-0.000139*T)*T)+
                ((0.30188-0.000344*T)+0.017998*dT)*dT )*dT/SAT_Const::Arcs;
    double z     =  zeta + ( (0.79280+0.000411*T)+0.000205*dT)*dT*dT/SAT_Const::Arcs;
    double theta =  ( (2004.3109-(0.85330+0.000217*T)*T)-
                ((0.42665+0.000217*T)+0.041833*dT)*dT )*dT/SAT_Const::Arcs;

    // Precession matrix
    return (r_z(-z) * r_y(theta) * r_z(-zeta));

 }