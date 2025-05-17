// $Header$
/**
 * @file accel.h
 * @author Marcos Gómez Robres
 * @brief Esta es la cabecera de la clase accel, esta calcula la aceleración de la tierra orbitando un satélite.
 * @version 0.1
 * @date 2025-04-30
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#ifndef _ACCEL_
#define _ACCEL_

#include "matrix.h"
#include "iers.h"
#include "timediff.h"
#include "global.h"
#include "precMatrix.h"
#include "nutMatrix.h"
#include "poleMatrix.h"
#include "JPL_Eph.h"
#include "accelPointMass.h"
#include "accelHarmonic.h"
#include "sat_const.h"
#include "mjday_tdb.h"
#include "ghaMatrix.h"
/**
 * @brief 
 * 
 * @param x Fecha terrestre
 * @param y Estado del satéliteq
 * @return Matrix aceleración del satélite
 */
Matrix accel(double x, Matrix y);

#endif