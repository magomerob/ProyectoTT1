// $Header$
/**
 * @file varEqn.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase poleMatrix, esta calcula la Transformación del sistema del meridiano de Greenwich a coordenadas tangentes locales
 * @date 2025-05-14
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#ifndef _VAREQN_
#define _VAREQN_
#include "matrix.h"
#include "iers.h"
#include "timediff.h"
#include "precMatrix.h"
#include "nutMatrix.h"
#include "poleMatrix.h"
#include "g_accelHarmonic.h"
#include "global.h"
#include "ghaMatrix.h"
#include "sat_const.h"

/**
 * @brief Esta es la cabecera de la clase varEqn, esta calcula la ecuación de variación
 * 
 * @param x tiempo desde la época
 * @param yPhi (6+36)-dim vector con el estado.
 * @return Matrix derivada de yPhi
 */

Matrix varEqn(double x, Matrix yPhi);

#endif