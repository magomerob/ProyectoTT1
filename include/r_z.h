// $Header$
/**
 * @file r_z.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase r_z, esta implementa un el cálculo de la matriz de rotación en el eje z para un ángulo cualquiera.
 * @date 2025-04-14
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#ifndef _rz_
#define _rz_
#include <cmath>
#include "matrix.h"
using namespace std;

/**
 * @brief Calcula la matriz de rotación en el eje z para un ángulo cualquiera.
 * 
 * @param angle Ángulo de rotación.
 * @return Matrix Matriz de rotación.
 */
Matrix r_z(double angle);

#endif
