// $Header$
/**
 * @file eqnEquinox.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase poleMatrix, esta calcula la ecuación de los equinoccios.
 * @date 2025-05-13
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#ifndef _EQNEQUINOX_
#define _EQNEQUINOX_
#include "nutAngles.h"
#include "meanObliquity.h"
#include <cmath>
/**
 * @brief calcula la ecuación de los equinoccios
 * 
 * @param Mjd_TT fecha juliana modificada
 * @return double ecuación de los equinoccios
 */

double eqnEquinox(double Mjd_TT);

#endif
