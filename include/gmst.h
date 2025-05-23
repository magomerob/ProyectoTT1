// $Header$
/**
 * @file gmst.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase gmst, esta implementa el tiempo medio Greenwich Sidereal.
 * @date 2025-05-07
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

 #ifndef _GMST_
 #define _GMST_
 #define _USE_MATH_DEFINES
 #include <cmath>
 #include "frac.h" 
 using namespace std;
 
 /**
  * @brief Calcula el tiempo medio Greenwich Sidereal.
  * 
  * @param Mjd_UT1 Fecha juliana modificada UT1.
  * @return double GMST en radianes.
  */
 double gmst(double Mjd_UT1);
 
 #endif
 