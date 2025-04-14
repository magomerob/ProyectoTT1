// $Source$
/**
 * @file position.cpp
 * @author Marcos Gómez Robres
 * @brief Implementación del cálculo de la posición geodésica.
 * @version 0.1
 * @date 2025-04-14
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#include "../include/position.h"
#include "../include/sat_const.h"

Matrix Position(double lon, double lat, double h)
{
    double R_equ = SAT_Const().R_Earth;
    double f = SAT_Const().f_Earth;

    double e2 = f*(2.0-f);   // Square of eccentricity
    double CosLat = cos(lat);// Cosine of geodetic latitude
    double SinLat = sin(lat);

    //Vector posición
    Matrix r(3);
    double N = R_equ / sqrt(1.0-e2*SinLat*SinLat);

    r(1) =  (         N+h)*CosLat*cos(lon);
    r(2) =  (         N+h)*CosLat*sin(lon);
    r(3) =  ((1.0-e2)*N+h)*SinLat;

    return r;
    }