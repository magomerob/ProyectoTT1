// $Source$
/**
 * @file meanObliquity.cpp
 * @author Marcos Gómez Robres
 * @brief Implementación el calculo de la obliquidad media.
 * @version 0.1
 * @date 2025-04-14
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#include "../include/meanObliquity.h"
#include "../include/sat_const.h"

double MOblq(double Mjd_TT)
{
    const double MJD_J2000 = SAT_Const::MJD_J2000;
    const double Rad = SAT_Const::Rad;

    double T = (Mjd_TT-MJD_J2000)/36525;

    return Rad *( 84381.448/3600-(46.8150+(0.00059-0.001813*T)*T)*T/3600 );
}