// $Header$
/**
 * @file meanObliquity.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase meanObliquity, esta implementa un el cálculo de la obliquidad media de una órbita elíptica.
 * @date 2025-04-14
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#ifndef _meanObliquity_
#define _meanObliquity_
#include <cmath>

using namespace std;

/**
 * @brief Calcula la obliquidad media de una órbita elíptica.
 * 
 * @param Mjd_TT Fecha Juliana modificada (Tiempo terrestre).
 * @return double Obliquidad media.
 */
double MOblq(double Mjd_TT);

#endif
