// $Header$
/**
 * @file eccAnom.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase eccAnom, esta implementa un el cálculo de la anomalía excéntrica para órbitas elípticas.
 * @version 0.1
 * @date 2025-04-14
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#ifndef _eccAnom_
#define _eccAnom_
#define _USE_MATH_DEFINES
#include <cmath>
#include "matrix.h"

using namespace std;

/**
 * @brief Calcula de la anomalía excéntrica de una órbita elíptica.
 * 
 * @param M Anomalía media en radianes.
 * @param e Excentricidad de la órbita de 0 a 1.
 * @return double Anomalía excentrica en radianes.
 */
double eccAnom(double M, double e);

#endif
