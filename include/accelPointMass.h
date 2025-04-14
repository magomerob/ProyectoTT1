// $Header$
/**
 * @file accelPointMass.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase accelPointMass, esta implementa este objeto en C++.
 * @version 0.1
 * @date 2025-04-11
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#ifndef _AccelPointMass_
#define _AccelPointMass_

#include <cmath>
#include "matrix.h"

/**
 * @brief Calcula la aceleración perturbacional dado un punto de masa.
 * 
 * @param r Vector posición del satélite.
 * @param s Vector posición del punto de masa.
 * @param GM Coeficiente gravitacional del punto de masa.
 * @return Matrix& Vector aceleración.
 */
Matrix AccelPointMass (Matrix &r, Matrix &s, double GM);	
#endif
