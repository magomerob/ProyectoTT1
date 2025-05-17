// $Header$
/**
 * @file timediff.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase timediff, esta implementa una lista de diferencias horarias.
 * @date 2025-04-14
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#ifndef _timediff_
#define _timediff_

#include <tuple>

using namespace std;

/**
 * @brief Función que devuelve las diferencias horarias.
 * 
 * @param UT1_UTC 
 * @param TAI_UTC 
 * @return tuple<double, double, double, double, double> Diccionario con los resultados.
 */
tuple<double, double, double, double, double> timediff(double UT1_UTC, double TAI_UTC);
#endif
