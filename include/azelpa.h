// $Header$
/**
 * @file azelpa.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la función azelpa, esta calcula el azimuth, la elevación y los parciales a partir de las coordenadas tangenciales.
 * @version 0.1
 * @date 2025-04-14
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#ifndef _azelpa_
#define _azelpa_

#define USING_MATH_DEFINES

#include <cmath>
#include "matrix.h"
#include <tuple>

using namespace std;

/**
 * @brief Calcula el azimuth, la elevación y los parciales a partir de las coordenadas tangenciales.
 * 
 * @param s Vector posición del satélite.
 * @return tuple<double, double, Matrix&, Matrix&> Azimuth, Elevación, Parciales en x e y.
 */
tuple<double, double, Matrix, Matrix> azelpa(Matrix& s);

#endif