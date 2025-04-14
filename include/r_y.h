// $Header$
/**
 * @file r_y.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase r_y, esta implementa un el cálculo de la matriz de rotación en el eje y para un ángulo cualquiera.
 * @date 2025-04-14
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#ifndef _ry_
#define _ry_
#include <cmath>
#include "matrix.h"
using namespace std;

/**
 * @brief Calcula la matriz de rotación en el eje y para un ángulo cualquiera.
 * 
 * @param angle Ángulo de rotación
 * @return Matrix Matriz de rotación.
 */
Matrix r_y(double angle);

#endif
