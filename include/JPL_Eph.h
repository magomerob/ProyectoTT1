// $Header$
/**
 * @file JPL_Eph.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase gmst, esta calcula la posición ecuatorial del Sol, la Luna y nueve planetas principales utilizando las Efemérides del JPL.
 * @date 2025-05-07
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

 #ifndef _JPL_
 #define _JPL_
 #include <cmath>
 #include <tuple>
 #include "matrix.h"
 #include "global.h"
 #include "cheb3d.h"
 using namespace std;
 
 /**
  * @brief Calcula la posición ecuatorial del Sol, la Luna y nueve planetas principales utilizando las Efemérides del JPL.
  * @param Mjd_TBD Fecha juliana modificada UT1.
  * @return <Matrix,Matrix,Matrix,Matrix,Matrix,Matrix,Matrix,Matrix,Matrix,Matrix,Matrix> posiciones de los diferentes planetas.
  */
 tuple<Matrix,Matrix,Matrix,Matrix,Matrix,Matrix,Matrix,Matrix,Matrix,Matrix,Matrix> JPL_Eph(double Mjd_TDB);
 
 #endif
 