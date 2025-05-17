// $Header$
/**
 * @file measUpdate.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase measUpdate
 * @version 0.1
 * @date 2025-04-30
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#ifndef _MEASUPDATE_
#define _MEASUPDATE_
#include "matrix.h"
#include <tuple>

using namespace std;

/**
 * @brief aplica la función meas update
 * 
 * @param x 
 * @param z 
 * @param g 
 * @param s 
 * @param G 
 * @param P 
 * @param n 
 * @return tuple<Matrix&, Matrix&, Matrix&> 
 */
tuple<Matrix, Matrix, Matrix> measUpdate(Matrix x, double z, double g, double s, Matrix G, Matrix P, int n);

#endif

