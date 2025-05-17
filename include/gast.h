// $Header$
/**
 * @file gast.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase gast, esta calcula el Tiempo sideral aparente de Greenwich
 * @date 2025-05-14
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#ifndef _GAST_
#define _GAST_

#include "gmst.h"
#include "eqnEquinox.h"
#include <cmath>

/**
 * @brief Esta es la cabecera de la clase gast, esta calcula el Tiempo sideral aparente de Greenwich
 * 
 * @param Mjd_UT1 tiempo en días julianos
 * @return double tiempo sideral aparente de Greenwich
 */
double gast(double Mjd_UT1);

#endif