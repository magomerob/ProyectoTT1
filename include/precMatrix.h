// $Header$
/**
 * @file precMatrix.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase precMatrix, esta implementa la transformación de la matriz de precesión.
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

 #ifndef _PRECMATRIX_
 #define _PRECMATRIX_
 #include <cmath>
 #include "matrix.h"
 #include "sat_const.h"
 #include "r_x.h"
 #include "r_y.h"
 #include "r_z.h"
 using namespace std;
 
 /**
  * @brief Calcula la transformación de la matriz de precesión.
  * 
  * @param Mjd_1 Fecha juliana modificada.
  * @param Mjd_2 Fecha juliana modificada.
  * @return Matrix matriz de precesión.
  */
 Matrix precMatrix(double Mjd_1, double Mjd_2);
 
 #endif
 