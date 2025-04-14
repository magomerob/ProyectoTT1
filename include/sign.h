// $Header$
/**
 * @file sign.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase sign, esta implementa una función que devuelve el valor absoluto del primer número con el signo de otro.
 * @date 2025-04-14
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#ifndef _sign_
#define _sign_
#include <cmath>

using namespace std;

/**
 * @brief Devuelve el valor absoluto del primer número con el signo del otro.
 * 
 * @param a Número del que se conserva el valor absoluto.
 * @param b Número del que se conserva el signo
 * @return double Valor absoluto de a con el signo de b.
 */
double sign(double a, double b);

#endif

