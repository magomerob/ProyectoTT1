// $Source$
/**
 * @file mjday_tdb.cpp
 * @author Marcos Gómez Robres
 * @brief Implementación del cálculo de la fecha juliana para tiempo dinámico baricentrico.
 * @version 0.1
 * @date 2025-04-14
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#include "../include/mjday_tdb.h"

double Mjd_TDB(double Mjd_TT)
{

    double T_TT = (Mjd_TT - 51544.5)/36525;

    double ret = Mjd_TT + ( 0.001658*sin(628.3076*T_TT + 6.2401)
                 +   0.000022*sin(575.3385*T_TT+4.2970)
                 +   0.000014*sin(1256.6152*T_TT + 6.1969)
                 +   0.000005*sin(606.9777*T_TT+4.0212)
                 +   0.000005*sin(52.9691*T_TT+0.4444)
                 +   0.000002*sin(21.3299*T_TT+5.543)
                 +   0.000010*sin(628.3076*T_TT+4.2490) )/86400;

    return ret;
}