// $Header$
/**
 * @file cheb3d.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase cheb3d, esta implementa una aproximación de chebyshev.
 * @version 0.1
 * @date 2025-04-11
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#ifndef _cheb3d_
#define _cheb3d_

#include <cmath>
#include "matrix.h"

/**
 * @brief Calcula la aproximación de chebyshev de vectores tridimensionales.
 * 
 * @param t Tiempo.
 * @param N Número de coeficientes.
 * @param Ta Inicio intervalo.
 * @param Tb Fin intervalo.
 * @param Cx Vector con coeficiente en coordenada x.
 * @param Cy Vector con coeficiente en coordenada y.
 * @param Cz Vector con coeficiente en coordenada z.
 * @return Matrix Vector con la aproximación.
 */
Matrix Cheb3D (double t,int N, double Ta, double Tb, Matrix &Cx, Matrix &Cy, Matrix &Cz);	
#endif
