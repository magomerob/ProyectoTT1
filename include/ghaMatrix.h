// $Header$
/**
 * @file ghaMatrix.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase ghaMatrix, esta calcula la transformación del ecuador verdadero y equinoccio al ecuador terrestre y sistema del meridiano de Greenwich
 * @date 2025-05-14
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#ifndef _GHAMATRIX_
#define _GHAMATRIX_

#include "matrix.h"
#include "r_z.h"
#include "gast.h"

/**
 * @brief Esta es la cabecera de la clase ghaMatrix, esta calcula la transformación del ecuador verdadero y equinoccio al ecuador terrestre y sistema del meridiano de Greenwich
 * 
 * @param Mjd_UT1 tiempo en días julianos
 * @return Matrix matriz angular de Greenwich
 */
Matrix GHAMatrix(double Mjd_UT1);

#endif

