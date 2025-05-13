// $Header$
/**
 * @file poleMatrix.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase poleMatrix, esta calcula la transformación de coordenadas pseudofijas de la Tierra a coordenadas fijas de la Tierra para una fecha determinada
 * @date 2025-05-13
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#ifndef _POLE_MATRIX_
#define _POLE_MATRIX_
#include "matrix.h"
#include "r_x.h"
#include "r_y.h"

/**
 * @brief calcula la transformación de coordenadas pseudofijas de la Tierra a coordenadas fijas de la Tierra para una fecha determinada.
 * 
 * @param xp desplazamiento en el eje x
 * @param yp desplazamiento en el eje y
 * @return Matrix matriz de transformación
 */

Matrix poleMatrix(double xp,double yp);

#endif