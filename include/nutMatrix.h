// $Header$
/**
 * @file nutMatrix.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase poleMatrix, esta calcula la transformación del ecuador medio al verdadero y equinoccio
 * @date 2025-05-13
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#ifndef _NUT_MATRIX_
#define _NUT_MATRIX_
#include "matrix.h"
#include "meanObliquity.h"
#include "nutAngles.h"
#include "r_x.h"
#include "r_y.h"
#include "r_z.h"

/**
 * @brief calcula la transformación del ecuador medio al verdadero y equinoccio
 * 
 * @param Mjd_TT fecha juliana modificada.
 * @return Matrix matriz de transformación.
 */

Matrix nutMatrix(double Mjd_TT);
#endif