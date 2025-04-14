// $Header$
/**
 * @file r_x.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase r_x, esta implementa un el cálculo de la matriz de rotación en el eje x para un ángulo cualquiera.
 * @date 2025-04-14
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#ifndef _rx_
#define _rx_
#include <cmath>
#include "matrix.h"
using namespace std;

/**
 * @brief Calcula la matriz de rotación en el eje x para un ángulo cualquiera.
 * 
 * @param angle Ángulo de rotación
 * @return Matrix Matriz de rotación.
 */
Matrix r_x(double angle);

#endif
