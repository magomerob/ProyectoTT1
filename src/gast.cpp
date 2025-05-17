// $Source$
/**
 * @file gast.cpp
 * @author Marcos Gómez Robres
 * @brief Implementación de gast.
 * @version 0.1
 * @date 2025-05-14
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#include "../include/gast.h"
#define _USE_MATH_DEFINES

double gast(double Mjd_UT1)
{
    return fmod(gmst(Mjd_UT1) + eqnEquinox(Mjd_UT1),  M_PI*2);
}