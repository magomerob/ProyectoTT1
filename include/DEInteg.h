// $Header$
/**
 * @file DEInteg.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase DEInteg, esta aplica métodos de integración numérica
 * @version 0.1
 * @date 2025-05-17
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#ifndef _DEINTEG_
#define _DEINTEG_

#include "matrix.h"
#include "sat_const.h"
#include "sign.h"
#include <cmath>

/**
 * @brief aplica métodos de integración numérica.
 * 
 * @param func 
 * @param t 
 * @param tout 
 * @param relerr 
 * @param abserr 
 * @param n_eqn 
 * @param y 
 * @return Matrix 
 */
Matrix DEInteg(Matrix func(double t, Matrix y), double t, double tout, double relerr, double abserr, int n_eqn, Matrix& y);

#endif
