// $Header$
/**
 * @file nutAngles.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la función NutAngles, esta calcula la nutación.
 * @version 0.1
 * @date 2025-04-24
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#ifndef _NUTANGLES_
#define _NUTANGLES_
#include <cmath>
#include "matrix.h"
#include "SAT_Const.h"
#include <tuple>

using namespace std;

/**
 * @brief Calcula la nutación.
 * 
 * @param Mjd_TT Fecha Juliana en tiempo Terrestre.
 * @return tuple nutación en longitud y oblicuidad.
 */

tuple<double,double> NutAngles (double Mjd_TT);

#endif