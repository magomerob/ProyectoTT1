// $Header$
/**
 * @file accelHarmonic.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase accelHarmonic, esta implementa este objeto en C++.
 * @version 0.1
 * @date 2025-04-30
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

 #ifndef _ACCELHARMONIC_
 #define _ACCELHARMONIC_
 
 #include <cmath>
 #include "matrix.h"
 #include "legendre.h"

 /**
  * @brief Calcula la aceleración perturbacional dado un punto de masa.
  * 
  * @param r Vector posición del satélite.
  * @param s Vector posición del punto de masa.
  * @param n_max Grados máximos.
  * @param m_max Orden máximo.
  * @return Matrix Vector aceleración.
  */
 Matrix accelHarmonic (Matrix &r, Matrix &E, double n_max, double m_max);	
 #endif
 