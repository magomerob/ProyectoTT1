// $Header$
/**
 * @file iers.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la función legendre, esta calcula dichos polinomios.
 * @version 0.1
 * @date 2025-04-15
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#ifndef _legendre_
#define _legendre_
#include <cmath>
#include "matrix.h"
#include <tuple>

using namespace std;

/**
 * @brief Calcula los polinomios de Legendre.
 * 
 * @param n Grado del polinomio.
 * @param m Orden del polinomio.
 * @param fi Ángulo en radianes.
 * @return tuple Valor del polinomio de Legendre.
 */

tuple<Matrix, Matrix> legendre(int n, int m, double fi);

#endif