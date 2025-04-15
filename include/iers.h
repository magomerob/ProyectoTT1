// $Header$
/**
 * @file iers.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la función iers, esta maneja el tiempo iers y el moimiento polar.
 * @version 0.1
 * @date 2025-04-15
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#ifndef _iers_
#define _iers_

#define USING_MATH_DEFINES

#include <cmath>
#include "matrix.h"
#include <tuple>
#include "sat_const.h"

using namespace std;

/**
 * @brief Calcula el tiempo iers y el movimiento polar.
 * 
 * @param eop Matriz de datos EOP.
 * @param Mjd_UTC Tiempo en MJD UTC.
 * @param interp Interpolación a realizar.
 * @return tuple<double, double, double, double, double, double, double, double, double> Datos iers.
 */
tuple<double, double, double, double, double, double, double, double, double> iers(Matrix& eop, double Mjd_UTC, char interp);

#endif